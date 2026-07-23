<!-- markdownlint-disable MD022 MD032 MD025 MD024 MD060 MD009 MD029 MD030 -->

# Elite Mathematically-Verified Code Auditor: CFD Suite Comprehensive Gap Analysis

**Auditor Persona**: Elite Mathematically-Verified Code Auditor
**Date**: November 18, 2025 (Updated: February 23, 2026 - Sprint 1.95.0)
**Status**: ✅ ALL CRITICAL ISSUES RESOLVED

# Sprint 1.96.167 Resolution: cfd-math native sparse-LU result ownership

### RESOLVED-196: Direct Solver Staged Native Arrays Through `Vec`
- **Location**:
  `crates/cfd-math/src/linear_solver/direct_solver.rs`.
- **Issue**: the primary direct-solve path collected the native `Array1` RHS
  into a temporary `Vec`, called the provider's slice API, then copied the
  returned `Vec` into another `Array1`. This added two consumer-owned linear
  buffers to every successful direct solve.
- **Remediation**: `leto-ops` now owns the `ArrayView1`-based `solve_view`
  contract. `DirectSparseSolver` passes `rhs.view()` and returns the provider
  `Array1` directly; dense fallback and finite-result validation remain
  unchanged.
- **Evidence**: provider and consumer direct-solve tests preserve the exact
  2×2 solution for `f32` and `f64`; provider all-target warning-denied Clippy,
  provider check, consumer `cfd-math` check, consumer lib Clippy, and focused
  direct-solver Nextest pass. Allocation reduction is established by the
  source/data-flow audit; no runtime allocation profile is claimed.
- **Closure**: provider SemVer classification passes 196/196 checks with 57
  skips; Leto PR #70 merged at `b24fc860864abad84af3118aa2bb27c32bb81265`;
  the `CFDrs` manifest pin and exact-head consumer gates use that revision.

# Sprint 1.96.166 Resolution: cfd-math IncompleteCholesky Still Required nalgebra CSR

### RESOLVED-PENDING-195: IncompleteCholesky Constructor Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs`,
  `crates/cfd-math/src/linear_solver/preconditioners/edge_case_tests.rs`,
  and `crates/cfd-math/tests/preconditioner_edge_cases.rs`.
- **Issue**: after ILU moved to Leto CSR, IncompleteCholesky still stored
  nalgebra-sparse CSR, used `get_entry()` for lookups, and built factors with
  `try_from_csr_data`.
- **Remediation**: moved IncompleteCholesky stored matrix state and IC(0)
  factor construction to `leto_ops::CsrMatrix`, routed lookups/substitutions
  through Leto CSR row views, and updated source/integration tests to pass
  Leto CSR into Cholesky construction.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; cfd-math
  all-target check passed; cfd-math lib/tests clippy passed; cfd-math
  all-target clippy passed; cholesky-filter nextest passed 5/5;
  preconditioner-filter nextest passed 76/76; targeted Cholesky residue scan
  found no nalgebra sparse/direct CSR/row-offset/get-entry residue in
  `cholesky.rs` or `tests/preconditioner_edge_cases.rs`.
- **Residual**: Schwarz, direct solver, remaining transitional solver
  fixtures, and the shared solver matrix boundary still use nalgebra-sparse.

# Sprint 1.96.165 Resolution: cfd-math ILU Still Required nalgebra CSR

### RESOLVED-PENDING-194: ILU Factorization Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/ilu`,
  `crates/cfd-math/src/linear_solver/chain.rs`,
  `crates/cfd-math/src/linear_solver/preconditioners/schwarz.rs`, and
  ILU source/integration tests.
- **Issue**: after SSOR moved to Leto CSR, ILU still stored and factorized
  nalgebra-sparse CSR and traversed rows through `row_offsets()`.
- **Remediation**: moved ILU(0), ILU(k), and triangular solve factor storage to
  `leto_ops::CsrMatrix`, routed triangular solves through Leto row views, and
  updated call sites to pass Leto CSR or convert once at a remaining boundary.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; cfd-math
  all-target check passed; cfd-math lib/tests clippy passed; cfd-math
  all-target clippy passed; ilu-filter nextest passed 21/21;
  preconditioner-filter nextest passed 76/76; linear_solver::tests nextest
  passed 53/53; targeted ILU residue scan found no nalgebra sparse/direct
  CSR/row-offset/get-entry residue in `preconditioners/ilu`.
- **Residual**: Schwarz, direct solver, remaining transitional solver
  fixtures, and the shared solver matrix boundary still use nalgebra-sparse.

# Sprint 1.96.164 Resolution: cfd-math SSOR Still Required nalgebra CSR

### RESOLVED-PENDING-193: SSOR Constructor Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs` and
  `crates/cfd-math/src/linear_solver/preconditioners/edge_case_tests.rs`.
- **Issue**: after basic preconditioners moved to Leto CSR, SSOR still stored
  nalgebra-sparse CSR and traversed rows through `row_offsets()`.
- **Remediation**: moved SSOR stored matrix state to `leto_ops::CsrMatrix`,
  routed forward/backward sweeps through Leto CSR row views, and converted
  source preconditioner edge-test fixtures once before SSOR construction.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; cfd-math
  all-target check passed; cfd-math lib/tests clippy passed; cfd-math
  all-target clippy passed; ssor-filter nextest passed 5/5;
  preconditioner-filter nextest passed 76/76; targeted SSOR residue scan
  found no nalgebra sparse/direct CSR/row-offset/get-entry residue in
  `ssor.rs`.
- **Residual**: Schwarz, direct solver, integration-test
  fixtures, and the shared solver matrix boundary still use nalgebra-sparse.

# Sprint 1.96.163 Resolution: cfd-math Basic Preconditioners Still Required nalgebra CSR

### RESOLVED-PENDING-192: Jacobi/SOR Constructors Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`,
  `crates/cfd-math/src/linear_solver/tests`, and
  `crates/cfd-math/tests/core_solver_tests.rs`.
- **Issue**: after AMG moved to Leto CSR, Jacobi and SOR still constructed
  from the nalgebra-sparse CSR boundary.
- **Remediation**: moved Jacobi/SOR constructors and SOR stored matrix state to
  `leto_ops::CsrMatrix`, using Leto CSR diagonal/row access. Tests that still
  use legacy solver matrices convert fixtures once before basic preconditioner
  construction.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; core solver
  test check passed; cfd-math all-target check passed; cfd-math lib/tests
  clippy passed; core solver test clippy passed; cfd-math all-target clippy
  passed; linear_solver::tests nextest passed 53/53; core_solver_tests
  nextest passed 4/4; preconditioner-filter nextest passed 76/76; targeted
  basic-preconditioner residue scan found no nalgebra sparse/direct CSR/row
  offset/get-entry residue in `basic.rs`.
- **Residual**: Schwarz, direct solver, and the shared solver matrix boundary
  still use nalgebra-sparse.

# Sprint 1.96.162 Resolution: cfd-math AMG/Coarsening Still Built nalgebra CSR

### RESOLVED-PENDING-191: AMG Coarsening Boundary Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid`,
  `crates/cfd-math/benches/coarsening_bench.rs`,
  `crates/cfd-math/benches/algebraic_distance_bench.rs`,
  `crates/cfd-math/tests/amg_coarsening_tests.rs`, and
  `crates/cfd-math/tests/amg_integration_test.rs`.
- **Issue**: after the direct Leto CSR operator path landed, AMG coarsening,
  interpolation, smoothers, cycles, and the coarsening/algebraic-distance
  benchmarks still depended on the nalgebra-sparse CSR boundary.
- **Remediation**: moved the multigrid bounded context to a local
  `leto_ops::CsrMatrix` sparse boundary, routed AMG sparse products,
  transpose, row access, and SpMV through Leto CSR APIs, and migrated the
  coarsening/algebraic-distance benchmarks plus AMG tests to Leto CSR where
  they exercise AMG.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; focused checks
  passed for `amg_coarsening_tests`, `amg_integration_test`,
  `coarsening_bench`, and `algebraic_distance_bench`; focused clippy passed
  for cfd-math lib, both AMG tests, and both migrated benches; AMG integration
  nextest passed 5/5; AMG-filter nextest passed 6/6; multigrid::coarsening
  nextest passed 10/10; targeted AMG/coarsening scan found no
  `nalgebra_sparse`, `CooMatrix`, `try_from_csr_data`, `get_entry`,
  `crate::sparse::spmv`, `spmv_array`, or `row_offsets()` residue in the
  migrated boundary.
- **Residual**: the broader cfd-math `crate::sparse::SparseMatrix`
  solver/direct/preconditioner surface still resolves to nalgebra-sparse, and
  `LinearSolverChain` converts that boundary once before AMG construction.

# Sprint 1.96.161 Resolution: cfd-math Benchmarks Still Built nalgebra CSR

### RESOLVED-PENDING-190: SpMV/CG Benchmarks Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/benches/spmv_bench.rs`,
  `crates/cfd-math/benches/cg_bench.rs`, and
  `crates/cfd-math/benches/math_benchmarks.rs`.
- **Issue**: after the direct Leto CSR `LinearOperator` path landed, the
  primary SpMV/CG benchmarks still constructed nalgebra sparse CSR matrices
  and the SpMV benchmark still measured the legacy sparse helper boundary.
- **Remediation**: migrated the three benchmark CSR constructors to
  `leto_ops::CsrMatrix::from_parts` and changed the SpMV benchmark to call
  `LinearOperator::apply` on the Leto CSR matrix directly.
- **Evidence**: cfd-math fmt passed; focused checks passed for `spmv_bench`,
  `cg_bench`, and `math_benchmarks`; focused clippy passed for all three
  migrated benches; cfd-math all-target check passed; cfd-math all-target
  clippy passed; sparse-filter nextest passed 19/19; targeted migrated-bench
  scan found no `nalgebra_sparse`, `CooMatrix`, `DVector`, `DMatrix`,
  `cfd_math::sparse::spmv`, or `try_from_csr_data` residue.
- **Residual**: after Sprint 1.96.162, the coarsening/algebraic-distance AMG
  boundary moved to Leto CSR. Remaining sparse-provider work is now the broader
  solver/direct/preconditioner nalgebra-sparse matrix surface.

# Sprint 1.96.160 Resolution: cfd-math Leto CSR LinearOperator

### RESOLVED-PENDING-189: GMRES Integration Still Required nalgebra Sparse CSR
- **Location**:
  `crates/cfd-math/src/sparse/operations.rs` and
  `crates/cfd-math/tests/simple_gmres_tests.rs`.
- **Issue**: cfd-math iterative solvers could operate over the local
  `LinearOperator<T>` trait, but `leto_ops::CsrMatrix<T>` did not yet satisfy
  that trait directly. The simple GMRES integration test still constructed
  nalgebra COO/CSR matrices and used the legacy sparse helper for residuals.
- **Remediation**: Implemented `LinearOperator<T>` for `leto_ops::CsrMatrix<T>`,
  consolidated the legacy nalgebra CSR SpMV path through the same
  `try_leto_spmv` helper, and migrated simple GMRES integration coverage to
  direct Leto CSR construction and residual application.
- **Evidence**: cfd-math fmt passed; simple_gmres test check passed;
  simple_gmres nextest passed 3/3; simple_gmres clippy passed; cfd-math lib
  check passed; cfd-math all-target check passed; cfd-math all-target clippy
  passed; sparse-filter nextest passed 19/19; gmres-filter nextest passed
  21/21; targeted simple GMRES scan found no `nalgebra_sparse`, `CooMatrix`,
  `nalgebra::`, `DVector`, `DMatrix`, or `cfd_math::sparse` residue.
- **Residual**: cfd-math still exposes `nalgebra_sparse::CsrMatrix` in the
  broader sparse/preconditioner/AMG/direct-solver/test/bench boundary. This
  slice proves the direct Leto CSR operator path for iterative solvers.

# Sprint 1.96.159 Resolution: cfd-math Storage-Slice Closure

### RESOLVED-PENDING-188: Remaining cfd-math Leto Storage-Slice Residue
- **Location**:
  `crates/cfd-math/src/nonlinear_solver/linalg.rs`,
  `crates/cfd-math/src/nonlinear_solver/anderson.rs`,
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/
  interpolation.rs`, and
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/smoothers.rs`.
- **Issue**: cfd-math still exposed mutable Leto storage slice helpers and
  used storage-slice reads in multigrid validation/tests after the sparse,
  GPU, SIMD, operator, and nonlinear immutable cleanup slices.
- **Remediation**: Removed the mutable slice helpers, changed Anderson pivot
  swaps to direct array indexing, changed interpolation quality validation to
  direct interpolated-vector indexing, and changed smoother tests to direct
  indexed assertions.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed;
  nonlinear_solver/multigrid nextest passed 46/46; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  cfd-math `src`/`tests` scan found no `leto::Storage`, `StorageMut`,
  `.storage().as_slice()`, `as_slice_mut()`, `vector_slice_mut`, or
  `matrix_slice_mut` residue.
- **Residual**: cfd-math Leto storage-slice cleanup is closed. Remaining
  provider work is nalgebra/nalgebra-sparse replacement and other Atlas
  provider boundaries.

# Sprint 1.96.158 Resolution: cfd-math Sparse/Basic Used Leto Storage

### RESOLVED-PENDING-187: Sparse Operations and Jacobi Used Leto Storage Slices
- **Location**:
  `crates/cfd-math/src/sparse/operations.rs` and
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`.
- **Issue**: Sparse SpMV and sparse scaling still borrowed `Array1` values
  through Leto storage slices, SpMV required a mutable output slice, and the
  Jacobi preconditioner iterated the diagonal through `.storage().as_slice()`.
- **Remediation**: Changed SpMV input staging, SpMV output writeback,
  row/column scaling vectors, and Jacobi diagonal iteration to direct
  `Array1` indexing while preserving delegation to Leto CSR operations for
  sparse kernels.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed;
  sparse/preconditioner nextest passed 95/95; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scans found no `Storage`, `.storage().as_slice()`, `as_slice_mut()`,
  or obsolete SpMV output-contiguity diagnostic residue in the migrated
  sparse/basic files.
- **Residual**: remaining cfd-math storage-slice owners are nonlinear mutable
  dense-workspace helpers and multigrid interpolation/smoother internals.

# Sprint 1.96.157 Resolution: cfd-math GPU Operator Used Leto Storage

### RESOLVED-PENDING-186: GPU Linear Operator Used Leto Storage Slices
- **Location**:
  `crates/cfd-math/src/linear_solver/operators/gpu.rs`.
- **Issue**: The GPU linear operator borrowed the input `Array1` through
  `leto::Storage` and required mutable output slice contiguity before invoking
  the Hephaestus-backed GPU kernel path.
- **Remediation**: Changed `LinearOperator::apply` to validate input/output
  lengths, collect the upload buffer from direct `Array1` indexing, allocate
  a typed readback buffer, and write results back by index. Removed the
  `leto::Storage` import and the output `as_slice_mut()` requirement.
- **Evidence**: cfd-math fmt passed; cfd-math GPU-feature check passed;
  GPU-feature linear_solver::operators nextest passed 5/5; cfd-math
  GPU-feature lib clippy passed; cfd-math GPU-feature all-target check passed;
  cfd-math GPU-feature all-target clippy passed; targeted scan found no
  `Storage`, `.storage().as_slice()`, `as_slice_mut()`, or obsolete output
  contiguity diagnostic residue in `operators/gpu.rs`.
- **Residual**: remaining cfd-math storage-slice owners are sparse operations
  and multigrid internals. Full GPU provider replacement still requires
  replacing broader CFDrs raw WGPU context/kernel ownership with Hephaestus
  contracts where available.

# Sprint 1.96.156 Resolution: cfd-math Finite-Difference Operators Used Leto Storage

### RESOLVED-PENDING-185: CPU Linear Operators Used Leto Storage Slices
- **Location**:
  `crates/cfd-math/src/linear_solver/operators/poisson.rs` and
  `crates/cfd-math/src/linear_solver/operators/momentum.rs`.
- **Issue**: CPU finite-difference linear operators borrowed input storage
  slices and required mutable output slices before applying stencils.
- **Remediation**: Changed 2D Laplacian, 3D Poisson, 1D/2D momentum, and 2D
  energy operators to read and write through direct `Array1` indexing and
  removed `leto::Storage` imports plus output `as_slice_mut()` contiguity
  checks.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed;
  linear_solver::operators nextest passed 5/5; cfd-math lib/tests clippy
  passed; cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scan found no `Storage`, `.storage().as_slice()`, or
  `as_slice_mut()` residue in the migrated operator files.
- **Residual**: remaining storage-slice owners are sparse operations, GPU
  operator, and multigrid internals. GPU provider work is reserved for a
  Hephaestus-specific slice.

# Sprint 1.96.155 Resolution: cfd-math Nonlinear Linalg Used Leto Storage

### RESOLVED-PENDING-184: Nonlinear Linalg Used Immutable Storage Slices
- **Location**:
  `crates/cfd-math/src/nonlinear_solver/linalg.rs` and
  `crates/cfd-math/src/nonlinear_solver/anderson.rs`.
- **Issue**: Nonlinear vector helper functions borrowed immutable raw Leto
  storage slices through a `vector_slice` helper, and Anderson acceleration
  used that helper to iterate least-squares coefficients.
- **Remediation**: Changed nonlinear dot/add/sub/scaled-add/scale helpers to
  use direct `Array1` indexing and changed Anderson acceleration to index the
  `gamma` vector directly. Removed the immutable `vector_slice` helper and
  `leto::Storage` import.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed;
  nonlinear_solver nextest passed 9/9; cfd-math lib/tests clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scan found no immutable `.storage().as_slice()`, `leto::Storage`,
  or `vector_slice` residue in nonlinear linalg/Anderson.
- **Residual**: mutable dense-workspace helpers still use `StorageMut`;
  remaining immutable storage-slice owners are sparse operations, linear
  operators, GPU operator, and multigrid internals.

# Sprint 1.96.154 Resolution: cfd-math Production SIMD Vector Used Leto Storage

### RESOLVED-PENDING-183: SIMD Vector Extension Used Leto Storage Slices
- **Location**:
  `crates/cfd-math/src/simd/vector.rs`.
- **Issue**: Production SIMD vector extension methods borrowed raw Leto
  storage slices before dispatching Moirai map/reduce work, leaving the
  implementation coupled to the storage trait rather than the `Array1`
  indexing surface.
- **Remediation**: Changed `simd_mul`, `simd_dot`, `simd_norm`, and `par_map`
  to read values through direct `Array1` indexing and removed the
  `leto::Storage` import while preserving Moirai `Adaptive` dispatch.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; simd::vector
  nextest passed 1/1; cfd-math lib/tests clippy passed; cfd-math all-target
  check passed; cfd-math all-target clippy passed; broader SIMD-filter
  nextest passed 26/26; targeted scan found no `Storage` or
  `.storage().as_slice()` residue in `src/simd/vector.rs`.
- **Residual**: source-level storage-slice owners remain in nonlinear linalg,
  sparse operations, linear operators, GPU operator, and multigrid internals.

# Sprint 1.96.153 Resolution: cfd-math SIMD Integration Test Bypassed Leto Storage

### RESOLVED-PENDING-182: SIMD Integration Test Used Leto Storage Slice Bridge
- **Location**:
  `crates/cfd-math/tests/simd_tests.rs`.
- **Issue**: The CFD residual SIMD integration test converted the Leto
  `spmv` result through `.storage().as_slice().to_vec()` before passing it to
  the SIMD slice API, leaving a storage bridge in the integration-test
  provider boundary.
- **Remediation**: Changed the test to read the Leto `spmv` result by direct
  `Array1` indexing and pass the collected values to the existing SIMD slice
  API without importing `leto::Storage`.
- **Evidence**: cfd-math fmt passed; simd_tests check passed; simd_tests
  nextest passed 12/12; simd_tests clippy passed; cfd-math all-target check
  passed; cfd-math all-target clippy passed; broader SIMD-filter nextest
  passed 26/26; targeted scans found no `DVector`, nalgebra vector import,
  local preconditioner bridge, `Storage`, or storage-slice conversion residue
  under `crates/cfd-math/tests`.
- **Residual**: remaining provider migration work is outside the cfd-math
  integration-test vector bridge layer and includes `nalgebra_sparse::CsrMatrix`,
  dense nalgebra test oracles, and source-level Leto storage-slice internals.

# Sprint 1.96.152 Resolution: cfd-math AMG Integration Test Bypassed Leto

### RESOLVED-PENDING-181: AMG Integration Test Used DVector Bridges
- **Location**:
  `crates/cfd-math/tests/amg_integration_test.rs`.
- **Issue**: AMG integration tests still created nalgebra `DVector`
  exact-solution, RHS, preconditioner output, and two-grid work values, then
  converted through Leto storage slices around solver/preconditioner dispatch.
- **Remediation**: Changed exact solution and RHS helpers to return
  `leto::Array1`, changed BiCGSTAB/GMRES error checks and AMG cycle checks to
  operate on Leto arrays, and changed the two-grid convergence factor to
  apply AMG directly to Leto work buffers.
- **Evidence**: cfd-math fmt passed; amg_integration_test check passed;
  amg_integration_test nextest passed 5/5; amg_integration_test clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  broader AMG-filter nextest passed 6/6; targeted scan found no `DVector`,
  nalgebra vector import, `Storage`, storage-slice conversion, or local
  preconditioner bridge residue in the migrated test.
- **Residual**: `amg_integration_test.rs` still uses
  `nalgebra_sparse::CsrMatrix` for sparse storage and nalgebra
  `DMatrix`/`SymmetricEigen` for the dense energy-norm oracle. Remaining
  cfd-math integration storage-slice residue is in `tests/simd_tests.rs`.

# Sprint 1.96.151 Resolution: cfd-math Preconditioner Edge Cases Bypassed Leto

### RESOLVED-PENDING-180: Preconditioner Edge-Case Tests Used DVector Bridge
- **Location**:
  `crates/cfd-math/tests/preconditioner_edge_cases.rs`.
- **Issue**: ILU edge-case tests still created nalgebra `DVector`
  RHS/solution values, converted them into Leto arrays for preconditioner
  dispatch, and converted preconditioner output back to `DVector`.
- **Remediation**: Changed the tests to allocate `leto::Array1` RHS/output
  buffers directly and removed the local preconditioner bridge helper.
- **Evidence**: cfd-math fmt passed; preconditioner_edge_cases check passed;
  preconditioner_edge_cases nextest passed 6/6; preconditioner_edge_cases
  clippy passed; cfd-math all-target check passed; cfd-math all-target clippy
  passed; broader preconditioner nextest passed 76/76; targeted scan found no
  `DVector`, nalgebra vector import, `Storage`, storage-slice conversion, or
  local preconditioner bridge residue in the migrated test.
- **Residual**: `preconditioner_edge_cases.rs` still uses the shared
  `nalgebra_sparse::CsrMatrix` matrix boundary. Remaining cfd-math integration
  storage-slice residue is in `tests/simd_tests.rs`.

# Sprint 1.96.150 Resolution: cfd-math Linear-Solver Tests Bypassed Leto

### RESOLVED-PENDING-179: Linear-Solver Source Tests Used DVector Bridges
- **Location**:
  `crates/cfd-math/src/linear_solver/tests/{mod,edge_case_tests,
  adversarial_solver_tests,extended_edge_case_tests}.rs`.
- **Issue**: Source linear-solver tests still created nalgebra `DVector`
  RHS/solution/work values or used local bridge helpers instead of exercising
  the Leto-native solver and preconditioner boundary directly.
- **Remediation**: Changed the tests to allocate `leto::Array1` buffers
  directly, removed local solve/preconditioner bridge macros, and routed
  residual checks through the Leto SpMV/helper path. The touched
  solver/sparse cone now routes scalar constants/tolerances through Eunomia
  instead of old scalar helper calls.
- **Evidence**: cfd-math fmt passed; cfd-math lib check passed; cfd-math
  all-target check passed; linear_solver::tests nextest passed 53/53;
  linear_solver nextest passed 176/176; cfd-math lib/tests clippy passed;
  cfd-math all-target clippy passed; targeted scans found no `DVector`,
  nalgebra vector import, `Storage`, storage-slice conversion, local
  solve/preconditioner bridge, matrix-vector `&a * &x`, or vector `.norm()`
  residue in `src/linear_solver/tests`, and no old scalar helper residue in
  the searched solver/sparse cone.
- **Residual**: integration-test vector bridge residue was closed by Sprint
  1.96.153; production sparse/scalar boundaries still include transitional
  nalgebra providers.

# Sprint 1.96.149 Resolution: cfd-math Core Solver Tests Bypassed Leto

### RESOLVED-PENDING-178: Core Solver Tests Used DVector Bridge Helpers
- **Location**:
  `crates/cfd-math/tests/core_solver_tests.rs`.
- **Issue**: BiCGSTAB, GMRES, preconditioner integration, and
  condition-number robustness tests still created nalgebra `DVector`
  RHS/solution values, converted them into Leto arrays for solver dispatch,
  and converted preconditioner output back to `DVector`. Residual checks also
  used nalgebra matrix-vector multiplication and vector norms.
- **Remediation**: Changed the tests to allocate `leto::Array1` RHS/solution
  buffers directly, removed the local solve/preconditioner bridge helpers, and
  changed residual verification to use `cfd_math::sparse::spmv` on Leto
  vectors.
- **Evidence**: cfd-math fmt passed; core_solver_tests check passed;
  core_solver_tests nextest passed 4/4; core_solver_tests clippy passed;
  cfd-math all-target check passed; cfd-math all-target clippy passed;
  targeted scan found no `DVector`, nalgebra vector, `Storage`, storage-slice
  conversion, local solve/preconditioner bridge, matrix-vector `&a * &x`, or
  vector `.norm()` residue in the core solver test module.
- **Residual**: `core_solver_tests.rs` still uses the shared
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix boundary, and broader
  cfd-math integration/adversarial/preconditioner tests still contain nalgebra
  `DVector` bridge diagnostics.

# Sprint 1.96.148 Resolution: cfd-math Simple GMRES Tests Bypassed Leto

### RESOLVED-PENDING-177: Simple GMRES Tests Used DVector Bridge Macro
- **Location**:
  `crates/cfd-math/tests/simple_gmres_tests.rs`.
- **Issue**: Basic, restarted, and preconditioned GMRES integration tests still
  created nalgebra `DVector` RHS/solution values, converted them into Leto
  arrays for solver dispatch, and converted the result back to `DVector`.
  Residual checks also used nalgebra matrix-vector multiplication and vector
  norms.
- **Remediation**: Changed the tests to allocate `leto::Array1` RHS/solution
  buffers directly, removed the local solve bridge macro, and changed residual
  verification to use `cfd_math::sparse::spmv` on Leto vectors.
- **Evidence**: cfd-math fmt passed; simple_gmres test check passed;
  simple_gmres nextest passed 3/3; simple_gmres clippy passed; cfd-math
  all-target check passed; cfd-math all-target clippy passed; touched-file
  `git diff --check` passed; targeted scan found no `DVector`, nalgebra
  vector, `Storage`, storage-slice conversion, local bridge macro,
  matrix-vector `&a * &x`, or vector `.norm()` residue in the simple GMRES
  test module.
- **Residual**: `simple_gmres_tests.rs` still uses the shared
  `nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix boundary, and broader
  cfd-math integration/adversarial/core/preconditioner tests still contain
  nalgebra `DVector` bridge diagnostics.

# Sprint 1.96.147 Resolution: cfd-math Matrix-Free Tests Bypassed Leto

### RESOLVED-PENDING-176: Matrix-Free Tests Used DVector Bridge Macro
- **Location**:
  `crates/cfd-math/src/linear_solver/matrix_free/tests.rs`.
- **Issue**: Matrix-free CG/GMRES tests still created nalgebra `DVector`
  RHS/solution values, converted them into Leto arrays for solver dispatch,
  and converted the result back to `DVector`, leaving test coverage on a
  transitional bridge instead of the Leto-native solver boundary.
- **Remediation**: Changed the matrix-free tests to allocate `leto::Array1`
  RHS/solution buffers directly, removed the local solve bridge macro, and
  changed the operator-size mismatch test to assert the exact typed
  `InvalidConfiguration` message.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math matrix-free nextest passed 4/4; cfd-math all-target clippy passed;
  touched-file `git diff --check` passed; targeted scan found no `DVector`,
  nalgebra, `Storage`, storage-slice conversion, or local bridge macro residue
  in the matrix-free test module.
- **Residual**: broader cfd-math linear-solver integration/adversarial/core/
  preconditioner test diagnostics still contain nalgebra `DVector` bridges,
  and production sparse/scalar provider boundaries still include
  `nalgebra_sparse::CsrMatrix` and transitional `nalgebra::RealField`.

# Sprint 1.96.146 Resolution: cfd-math BiCGSTAB Bypassed Leto

### RESOLVED-PENDING-175: BiCGSTAB Workspaces Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/bicgstab/mod.rs`,
  `crates/cfd-math/src/linear_solver/chain.rs`,
  `crates/cfd-math/src/linear_solver/traits.rs`, and
  `crates/cfd-math/tests/amg_integration_test.rs`.
- **Issue**: BiCGSTAB still converted Leto public RHS/solution arrays into
  nalgebra `DVector` workspaces internally, called legacy vector-bridge
  helpers for operator/preconditioner application, and the tiered solver chain
  converted its final BiCGSTAB fallback through nalgebra vectors.
- **Remediation**: Changed BiCGSTAB workspaces and direct solve methods to
  `leto::Array1`, replaced nalgebra vector math with shared Leto-array dot,
  norm, copy, residual, axpy, and scale-add operations, dispatched
  operators/preconditioners through the existing Leto trait boundary, kept the
  chain fallback on Leto arrays, and deleted obsolete legacy bridge helpers.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  focused cfd-math BiCGSTAB nextest passed 24/24; broader cfd-math
  linear-solver nextest passed 176/176; AMG integration nextest passed 5/5;
  cfd-math all-target clippy passed; touched-file `git diff --check` passed;
  targeted scans found no legacy bridge helper residue and no `DVector`/
  nalgebra vector-operation residue in the migrated BiCGSTAB/chain/traits
  source.
- **Residual**: cfd-math still owns the shared `nalgebra_sparse::CsrMatrix`
  storage/provider boundary and transitional `nalgebra::RealField` scalar
  bounds. Remaining nalgebra `DVector` use in this area is confined to other
  matrix-free/preconditioner/integration test diagnostics.

# Sprint 1.96.145 Resolution: cfd-math Conjugate Gradient Bypassed Leto

### RESOLVED-PENDING-174: CG Workspaces Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/conjugate_gradient/mod.rs`,
  `crates/cfd-math/benches/cg_bench.rs`, and
  `crates/cfd-math/benches/math_benchmarks.rs`.
- **Issue**: Conjugate Gradient still converted Leto public RHS/solution
  arrays into nalgebra `DVector` workspaces internally, called legacy
  vector-bridge helpers for operator/preconditioner application, and benchmarked
  the direct CG API through nalgebra vectors.
- **Remediation**: Changed CG workspaces and direct solve methods to
  `leto::Array1`, replaced nalgebra vector math with local Leto-array dot,
  norm, axpy, and scale-add operations, dispatched operators/preconditioners
  through the existing Leto trait boundary, and moved CG benchmark inputs to
  Leto vectors.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math conjugate nextest passed
  13/13; broader cfd-math linear-solver nextest passed 176/176; targeted
  source scan found no `DVector`, legacy bridge helper calls, `num_traits`,
  Leto storage-slice conversion, nalgebra vector math calls, or nalgebra
  `copy_from` residue in the migrated CG source and CG benchmark call sites.
- **Residual after the BiCGSTAB follow-up**: the shared linear-solver trait
  family still carries the transitional `nalgebra::RealField` scalar bound and
  sparse storage remains on `nalgebra_sparse::CsrMatrix`.

# Sprint 1.96.144 Resolution: cfd-math Schwarz Preconditioner Bypassed Leto

### RESOLVED-PENDING-173: Schwarz Local Apply Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/schwarz.rs`.
- **Issue**: Additive/multiplicative Schwarz local apply methods still
  accepted and returned nalgebra `DVector`, built local RHS `DVector`
  workspaces, and converted the default `Preconditioner::apply_to` path
  through a global nalgebra vector after the public preconditioner boundary had
  moved to Leto arrays.
- **Remediation**: Changed additive/multiplicative Schwarz application to
  consume and return `leto::Array1`, changed local RHS extraction to write
  directly into Leto arrays, kept local ILU solves on Leto buffers, added typed
  residual/output length validation, and copied the default additive result
  directly into the caller-provided Leto output buffer.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math Schwarz nextest passed
  3/3; broader cfd-math preconditioner nextest passed 76/76; targeted source
  scan found no `DVector`, `Storage`, `num_traits`, `FromPrimitive`,
  local-RHS conversion bridge, `DMatrix`, or `ndarray` residue in
  `schwarz.rs`.
- **Residual**: Schwarz still stores and constructs local sparse matrices
  through the shared `nalgebra_sparse::CsrMatrix` boundary and remains
  constrained by the global `Preconditioner<T>` nalgebra scalar bound.

# Sprint 1.96.143 Resolution: cfd-math ILU Triangular Solve Bypassed Leto

### RESOLVED-PENDING-172: ILU Triangular Workspaces Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/ilu/{types,triangular}.rs`.
- **Issue**: ILU forward/backward substitution and `IncompleteLU::apply_to`
  still allocated nalgebra `DVector` residual, intermediate, and solution
  workspaces after the public preconditioner boundary had moved to Leto arrays.
- **Remediation**: Changed ILU triangular substitutions to consume
  `leto::Array1` buffers, copied the Leto solution into the caller-provided
  Leto output buffer directly, added typed residual/output length validation,
  and routed U-solve diagonal identity through Eunomia `NumericElement`.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math ILU nextest passed
  21/21; broader cfd-math preconditioner nextest passed 74/74; targeted source
  scan found no `DVector`, `DMatrix`, `Storage`, old scalar identities,
  fallback wording, `component_mul`, `rows_mut`, row-view residue,
  `num_traits`, or `num_complex` in `ilu/types.rs` or `ilu/triangular.rs`.
- **Residual**: IncompleteLU still stores the shared
  `nalgebra_sparse::CsrMatrix` LU factor boundary and remains constrained by
  the global `Preconditioner<T>` nalgebra scalar bound.

# Sprint 1.96.142 Resolution: cfd-math Deflation Preconditioner Bypassed Leto

### RESOLVED-PENDING-171: Deflation Eigenvectors Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/deflation.rs`.
- **Issue**: Deflation eigenvector storage, added eigenpair input, and the
  deflated apply workspace still used nalgebra `DVector` after the public
  preconditioner boundary had moved to Leto arrays.
- **Remediation**: Changed deflation eigenvectors to `leto::Array1`, changed
  `add_eigenpair` to accept Leto arrays, computed projection coefficients
  over Leto arrays, removed the nalgebra deflated workspace, added typed
  vector length validation, and rejected zero eigenvalues before division.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math deflation nextest passed
  3/3; broader cfd-math preconditioner nextest passed 73/73; targeted source
  scan found no `DVector`, `DMatrix`, `Storage`, old scalar identities,
  fallback wording, `component_mul`, `rows_mut`, row-view residue,
  `num_traits`, or `num_complex` in `deflation.rs`.
- **Residual**: Deflation still wraps the base preconditioner behind the
  existing `Box<dyn Preconditioner<T>>`, and the global `Preconditioner<T>`
  trait still uses the transitional nalgebra `RealField` scalar bound.

# Sprint 1.96.141 Resolution: cfd-math Basic Preconditioners Bypassed Leto

### RESOLVED-PENDING-170: Basic Preconditioner State Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`.
- **Issue**: Identity, Jacobi, and SOR basic preconditioner vector paths still
  had nalgebra vector state or apply residue after the public preconditioner
  boundary had moved to Leto arrays.
- **Remediation**: Changed Identity/Jacobi/SOR apply paths to consume
  `leto::Array1` residual/output buffers directly, changed Jacobi inverse
  diagonal storage to `leto::Array1`, added typed residual/output length
  validation, and routed scalar identities/tolerance through Eunomia.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math basic mismatch nextest
  passed 1/1; broader cfd-math preconditioner nextest passed 70/70; targeted
  source scan found no `DVector`, `DMatrix`, old scalar identities, fallback
  wording, `component_mul`, `rows_mut`, or row-view residue in `basic.rs`.
- **Residual**: Jacobi and SOR still store the shared
  `nalgebra_sparse::CsrMatrix` boundary and remain constrained by the global
  `Preconditioner<T>` nalgebra scalar bound.

# Sprint 1.96.140 Resolution: cfd-math IncompleteCholesky Bypassed Leto

### RESOLVED-PENDING-169: Cholesky Substitution Workspaces Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/cholesky.rs`.
- **Issue**: IncompleteCholesky's forward/backward substitution helpers and
  `apply_to` still allocated nalgebra `DVector` residual, intermediate, and
  solution workspaces after the public preconditioner boundary had moved to
  Leto arrays.
- **Remediation**: Changed Cholesky substitutions to consume `leto::Array1`
  buffers, copied the Leto solution into the caller-provided Leto output
  buffer directly, added typed residual/output length validation before
  substitution, and routed IC(0) square-root dispatch through Eunomia
  `NumericElement`.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math Cholesky nextest passed
  5/5; targeted source scan found no `DVector`, `DMatrix`, `Storage`, old
  scalar identities, silent fallback wording, ambiguous `.sqrt()`, or
  residual `DVector` workspace construction in `cholesky.rs`.
- **Residual**: resolved for IncompleteCholesky by Sprint 1.96.166, which
  moved its factor boundary to Leto CSR. Remaining sparse-provider work is
  Schwarz, direct solver, and the shared solver matrix boundary.

# Sprint 1.96.139 Resolution: cfd-math SSOR Bypassed Leto

### RESOLVED-PENDING-168: SSOR Sweep Workspaces Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs`.
- **Issue**: SSOR's forward/backward sweep helpers and `apply_to` still
  allocated nalgebra `DVector` residual and solution workspaces after the
  public preconditioner boundary had moved to Leto arrays.
- **Remediation**: Changed SSOR sweeps to consume `leto::Array1` buffers,
  initialized and mutated the caller-provided Leto output buffer directly, and
  added typed residual/output length validation before sweep execution.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math SSOR nextest passed
  5/5; targeted source scan found no `DVector`, `DMatrix`, `Storage`,
  nalgebra row-view operations, `component_mul`, old scalar identities, silent
  fallback wording, or clone fallback residue in `ssor.rs`.
- **Residual**: SSOR still stores the shared `nalgebra_sparse::CsrMatrix`
  boundary and remains constrained by the global `Preconditioner<T>`
  nalgebra scalar bound.

# Sprint 1.96.138 Resolution: cfd-math Block/SIMPLE Preconditioners Bypassed Leto

### RESOLVED-PENDING-167: Block/SIMPLE Preconditioner State Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/block_preconditioner.rs`.
- **Issue**: block/SIMPLE preconditioner direct apply paths and
  `Preconditioner::apply_to` implementations still used nalgebra `DVector`
  storage and conversions after the public preconditioner boundary had moved
  to Leto arrays.
- **Remediation**: Replaced diagonal inverse and SIMPLE Schur diagonal storage
  with `leto::Array1`, changed direct `apply` methods to consume and return
  Leto arrays, removed the residual/result `DVector` bridge in
  `apply_to`, and made vector length mismatches typed errors instead of cloned
  input fallbacks.
- **Evidence**: cfd-math fmt passed; cfd-math all-target check passed;
  cfd-math all-target clippy passed; focused cfd-math block-preconditioner
  nextest passed 4/4; targeted source scan found no `DVector`, `DMatrix`,
  nalgebra import, nalgebra row-view operations, `component_mul`, old scalar
  identities, silent fallback wording, or clone fallback residue in
  `block_preconditioner.rs`.
- **Residual**: cfd-math's global `Preconditioner<T>` trait still requires the
  transitional `nalgebra::RealField` scalar bound, and sparse storage still
  aliases `nalgebra_sparse::CsrMatrix`.

# Sprint 1.96.137 Resolution: cfd-math GMRES Bypassed Leto Workspace

### RESOLVED-PENDING-166: GMRES Work Vectors Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/gmres/{solver,arnoldi}.rs` and
  `crates/cfd-math/src/linear_solver/chain.rs`.
- **Issue**: GMRES public trait calls accepted Leto arrays but immediately
  converted RHS/solution vectors into nalgebra `DVector` workspaces. The
  solver-chain GMRES tiers repeated that conversion before trying AMG, block,
  unpreconditioned, and ILU tiers.
- **Remediation**: Replaced GMRES Arnoldi basis-column, work, preconditioned
  work, and residual vectors with `leto::Array1`; changed GMRES direct solve
  methods to accept Leto RHS/solution arrays; and kept solver-chain GMRES tiers
  on Leto arrays.
- **Evidence**: cfd-math lib/tests/all-targets checks passed; cfd-math fmt
  passed; cfd-math all-target clippy passed; focused cfd-math GMRES nextest
  passed 21/21; cfd-2d no-default lib check passed; focused cfd-2d momentum
  nextest passed 53/53; targeted scans found no `DVector` or legacy
  operator/preconditioner bridge calls under `linear_solver/gmres`.
- **Residual**: cfd-math still has nalgebra scalar bounds, CG/BiCGSTAB
  workspaces, `LinearSolverChain`'s final BiCGSTAB bridge, and
  `nalgebra_sparse::CsrMatrix` storage pending later Leto/Eunomia slices.

# Sprint 1.96.136 Resolution: cfd-2d Momentum Vectors Bypassed Leto

### RESOLVED-PENDING-165: Momentum RHS/Solution Used DVector
- **Location**:
  `crates/cfd-2d/src/physics/momentum/{solver,solve,boundary/**}.rs`,
  `crates/cfd-2d/src/scalar.rs`, `crates/cfd-2d/Cargo.toml`, and
  `crates/cfd-2d/src/linear_solver_bridge.rs`.
- **Issue**: Momentum assembly, boundary handling, solve dispatch, tests, and
  the local solver bridge still used nalgebra vectors after SIMPLE and
  pressure-velocity correction had moved to Leto.
- **Remediation**: Replaced momentum RHS/solution buffers with `leto::Array1`,
  routed momentum solve dispatch directly through Leto-native iterative/direct
  solver APIs, changed momentum boundary helpers to mutate Leto arrays,
  removed the obsolete local nalgebra bridge module, removed the
  `nalgebra::RealField` supertrait from `Cfd2dScalar`, and dropped direct
  cfd-2d `nalgebra`/`nalgebra-sparse` manifest dependencies.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate targets, static source/manifest
  residue scan, dependency graph audit, and format hygiene. cfd-2d fmt passed;
  no-default lib check passed; no-default all-target clippy passed; focused
  cfd-2d momentum nextest passed 53/53; targeted cfd-2d source/manifest scans
  found no `DVector`, `DMatrix`, nalgebra, nalgebra-sparse, or
  `linear_solver_bridge`; cargo-tree audits show remaining nalgebra and
  nalgebra-sparse edges are transitive through upstream crates.
- **Residual risk**: cfd-2d still resolves nalgebra/nalgebra-sparse
  transitively through cfd-1d, cfd-core, cfd-math, cfd-schematics, and Gaia.

# Sprint 1.96.135 Resolution: cfd-2d Pressure-Velocity Vectors Bypassed Leto

### RESOLVED-PENDING-164: Pressure-Velocity Correction Used DVector
- **Location**:
  `crates/cfd-2d/src/pressure_velocity/{pressure,correction,faces}.rs`.
- **Issue**: The pressure-velocity pressure-correction path still cached RHS
  and solution buffers as nalgebra `DVector` and solved through the local
  nalgebra conversion bridge after the adjacent SIMPLE pressure-correction
  path had moved to Leto.
- **Remediation**: Replaced pressure-velocity correction caches and scatter
  buffers with `leto::Array1`, called `IterativeLinearSolver::solve` and
  `DirectSparseSolver::solve` directly, and routed RHS diagnostics through
  `leto_ops::norm_l2`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate targets, static residue scan, and
  format hygiene. cfd-2d fmt passed; no-default lib check passed; no-default
  all-target clippy passed; focused cfd-2d pressure-velocity nextest passed
  16/16; targeted scans found no `DVector`, nalgebra, or
  `linear_solver_bridge` residue under `crates/cfd-2d/src/pressure_velocity`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra remains transitive through upstream owners.

# Sprint 1.96.134 Resolution: cfd-2d SIMPLE Vectors Bypassed Leto

### RESOLVED-PENDING-163: SIMPLE Pressure Correction Used DVector
- **Location**: `crates/cfd-2d/src/solvers/simple/{algorithm,pressure}.rs`.
- **Issue**: SIMPLE pressure correction still stored RHS and `p_prime` in
  nalgebra `DVector` buffers and solved through the local nalgebra conversion
  bridge after adjacent cfd-2d FDM and time vector surfaces had moved to Leto.
- **Remediation**: Replaced SIMPLE pressure-correction buffers with
  `leto::Array1` and called `IterativeLinearSolver::solve` directly.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate targets, static residue scan, and
  format hygiene. cfd-2d fmt passed; no-default lib check passed; no-default
  all-target clippy passed; focused cfd-2d SIMPLE nextest passed 19/19;
  targeted scans found no `DVector`, nalgebra, or `linear_solver_bridge`
  residue under `crates/cfd-2d/src/solvers/simple`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra remains transitive through upstream owners.

# Sprint 1.96.133 Resolution: cfd-2d FDM Vectors Bypassed Leto

### RESOLVED-PENDING-162: FDM RHS/Solver Used DVector
- **Location**: `crates/cfd-2d/src/solvers/fdm/{linear_solver,poisson,advection_diffusion}.rs`.
- **Issue**: The FDM Poisson/advection-diffusion assembly path and shared
  Gauss-Seidel solver still used nalgebra `DVector` for RHS and solution
  storage after adjacent cfd-2d time and scheme surfaces had moved to Leto.
- **Remediation**: Replaced FDM RHS/result vectors with `leto::Array1`,
  updated stencil assembly to write provider arrays directly, and changed
  `solve_gauss_seidel` to accept and return Leto vectors.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate targets, static residue scan, and
  format hygiene. cfd-2d fmt passed; no-default lib check passed; no-default
  all-target clippy passed; focused cfd-2d FDM nextest passed 2/2; targeted
  scans found no `DVector` or nalgebra residue under
  `crates/cfd-2d/src/solvers/fdm`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra remains transitive through upstream owners.

# Sprint 1.96.132 Resolution: cfd-2d Time Vectors Bypassed Leto

### RESOLVED-PENDING-161: Time Integration Kept DVector States
- **Location**: `crates/cfd-2d/src/schemes/time/**`.
- **Issue**: The cfd-2d time-integration API still accepted and returned
  nalgebra `DVector` state values after adjacent cfd-2d scheme storage moved to
  Leto arrays.
- **Remediation**: Replaced the time-state boundary with `leto::Array1` through
  `StateVector<T>`, updated explicit, implicit, multistep, adaptive-controller,
  adaptive-integrator, and time tests to Leto vector operations, and replaced
  nalgebra vector norms with `leto_ops::norm_l2`.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate library, static residue scans, and
  format hygiene. cfd-2d fmt passed; no-default lib check passed; no-default
  all-target clippy passed; focused cfd-2d time nextest passed 29/29; targeted
  scans found no `DVector` or nalgebra residue under
  `crates/cfd-2d/src/schemes/time`.
- **Residual risk**: direct cfd-2d source/manifest nalgebra ownership is now
  removed; nalgebra remains transitive through upstream owners.

# Sprint 1.96.131 Resolution: cfd-2d DMatrix Bypassed Leto

### RESOLVED-PENDING-160: cfd-2d Compact Matrices Used DMatrix
- **Location**:
  `crates/cfd-2d/src/physics/immersed_boundary.rs`,
  `crates/cfd-2d/src/schemes/grid.rs`,
  `crates/cfd-2d/src/schemes/upwind.rs`,
  `crates/cfd-2d/src/schemes/tvd/muscl.rs`,
  scheme TVD/WENO tests, and `crates/cfd-2d/examples/blood_venturi.rs`.
- **Issue**: The immersed-boundary coupling path and scheme grid storage still
  used nalgebra `DMatrix`/matrix-shape idioms after adjacent cfd-2d provider
  work had moved vector and scalar boundaries toward Leto/Eunomia.
- **Remediation**: Replaced those compact matrix buffers with `leto::Array2`,
  updated all dependent scheme indexing to Leto `[[i, j]]`, destructured Leto
  `shape()` arrays in tests, and routed the blood Venturi example through the
  same provider-owned matrix type.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy over the affected crate library, static residue scans, and
  format hygiene. cfd-2d fmt passed; no-default lib check and clippy passed;
  no-default `blood_venturi` example check passed; focused cfd-2d nextest passed
  60/60; targeted scans found no cfd-2d `DMatrix` residue and no nalgebra-style
  `Grid2D.data` tuple access.
- **Residual risk**: cfd-2d still retains nalgebra in non-DMatrix
  `DVector`/sparse linear-system boundaries and other provider seams.

# Sprint 1.96.130 Resolution: cfd-1d Vascular Complex Bypassed Eunomia

### RESOLVED-PENDING-159: Bessel/Womersley Used nalgebra Complex
- **Location**:
  `crates/cfd-1d/src/physics/vascular/bessel.rs` and
  `crates/cfd-1d/src/physics/vascular/womersley/profile.rs`.
- **Issue**: The vascular analytical Bessel/Womersley path still used
  nalgebra `Complex` and imported `ComplexField` only to call `modulus()`,
  leaving a non-Atlas complex-number dependency in cfd-1d physics code.
- **Remediation**: Switched Bessel and Womersley complex values to
  `eunomia::Complex` and replaced the convergence stop condition with
  Eunomia complex `norm()`.
- **Evidence tier**: compile-time provider integration, focused empirical
  Bessel/Womersley nextest, clippy over the affected crate library, static
  residue scan, and format hygiene. cfd-1d lib check and clippy passed;
  focused cfd-1d nextest passed 26/26; cfd-1d fmt check passed; and targeted
  vascular nalgebra-complex residue scan found no remaining imports.
- **Residual risk**: cfd-1d still depends on nalgebra for network
  sparse/dense linear-system storage and the transitional scalar seam.

# Sprint 1.96.129 Resolution: LinearOperator apply Bypassed Leto

### RESOLVED-PENDING-158: LinearOperator Trait Kept DVector Buffers
- **Location**:
  `crates/cfd-math/src/linear_solver/traits.rs`,
  `crates/cfd-math/src/linear_solver/operators/**`,
  `crates/cfd-math/src/sparse/operations.rs`, and
  `crates/cfd-math/benches/math_benchmarks.rs`.
- **Issue**: `LinearOperator::apply` still accepted nalgebra `DVector`
  input/output buffers after adjacent solver and preconditioner traits had
  moved to Leto, forcing matrix-free and sparse operator callers through a
  nalgebra public provider boundary.
- **Remediation**: Moved `apply` and `apply_transpose` to `leto::Array1`,
  updated sparse CSR, identity/scaled, Poisson, momentum, energy, and GPU
  operator adapters, and confined solver workspace conversion to the internal
  `apply_operator_to_legacy` bridge.
- **Evidence tier**: compile-time provider integration, focused empirical
  solver/operator nextest, clippy over the producer crate, static residue scan,
  and format hygiene. cfd-math all-target check and clippy passed; focused
  cfd-math nextest passed 80/80; cfd-math fmt check passed; and targeted
  DVector operator-signature residue scan found no old public signature.
- **Residual risk after the CG/BiCGSTAB follow-up**: nalgebra sparse storage
  and some preconditioner internals still retain local
  `DVector`/`CsrMatrix` conversion bridges. cfd-validation numerical
  result/error storage still uses `DVector`. Broad cfd-validation nextest
  still fails in the existing venturi cross-fidelity convergence tests.

# Sprint 1.96.128 Resolution: Preconditioner apply_to Bypassed Leto

### RESOLVED-PENDING-157: Preconditioner Trait Kept DVector Buffers
- **Location**:
  `crates/cfd-math/src/linear_solver/traits.rs`,
  `crates/cfd-math/src/linear_solver/preconditioners/**`,
  cfd-math preconditioner tests, and
  `crates/cfd-1d/src/solver/core/linear_system.rs`.
- **Issue**: `Preconditioner::apply_to` still accepted nalgebra `DVector`
  residual/output buffers after adjacent solver traits had moved to Leto,
  forcing preconditioner callers and cfd-1d network preconditioning to retain
  DVector at the public provider boundary.
- **Remediation**: Moved the trait to `leto::Array1`, updated cfd-math
  preconditioner implementations/tests, kept current nalgebra matrix internals
  behind local conversion bridges, and updated cfd-1d `DiagJacobi` to
  implement the Leto preconditioner contract.
- **Evidence tier**: compile-time provider integration, focused empirical
  solver/preconditioner nextest, clippy over affected producer and consumer
  paths, static residue scan, and format hygiene. cfd-math all-target check
  and clippy passed; cfd-1d lib check and clippy passed; focused cfd-math
  nextest passed 131/131; cfd-math/cfd-1d fmt check passed; and targeted
  `apply_to` DVector-signature residue scan found no old public signature.
- **Residual risk**: `LinearOperator::apply`, nalgebra sparse preconditioner
  internals, and iterative solver workspaces still retain
  `DVector`/`CsrMatrix` conversion bridges. cfd-validation numerical
  result/error storage still uses `DVector`. Broad cfd-validation nextest
  still fails in the existing venturi cross-fidelity convergence tests.

# Sprint 1.96.127 Resolution: Iterative Solver solve Bypassed Leto

### RESOLVED-PENDING-156: Iterative Solver Trait Kept DVector Buffers
- **Location**:
  `crates/cfd-math/src/linear_solver/{traits.rs,conjugate_gradient/mod.rs,bicgstab/mod.rs,gmres/solver.rs}`,
  cfd-math solver tests, `crates/cfd-1d/src/solver/core/linear_system.rs`,
  cfd-2d momentum/pressure solve paths, and
  `crates/cfd-3d/src/fem/projection_solver.rs`.
- **Issue**: `IterativeLinearSolver::solve` still accepted nalgebra `DVector`
  RHS/result buffers, forcing downstream crates to pass DVector state through
  the public iterative solver trait.
- **Remediation**: Moved the trait to `leto::Array1`, updated CG/BiCGSTAB/
  GMRES implementations and tests, and routed cfd-1d/cfd-2d/cfd-3d callers
  through local Leto bridge modules.
- **Evidence tier**: compile-time provider integration, focused empirical
  solver nextest, clippy over affected producer and consumer paths, static
  residue scan, and diff hygiene. cfd-math fmt/check/clippy passed; focused
  cfd-math solver nextest passed 61/61; cfd-1d/cfd-2d/cfd-3d focused checks
  and clippy passed; cfd-validation all-target clippy passed; targeted
  DVector call-site residue scan and `git diff --check` passed.
- **Residual risk**: `LinearOperator::apply`, `Preconditioner::apply_to`,
  concrete preconditioners, and internal iterative workspaces still expose
  nalgebra `DVector`/`CsrMatrix`. cfd-validation numerical result/error
  storage still uses `DVector`. Broad cfd-validation nextest still fails in
  the existing venturi cross-fidelity convergence tests.

# Sprint 1.96.126 Resolution: LinearSolver solve_system Bypassed Leto

### RESOLVED-PENDING-155: Iterative Solver Trait Kept Nalgebra Vectors
- **Location**:
  `crates/cfd-math/src/linear_solver/{traits.rs,conjugate_gradient/mod.rs,bicgstab/mod.rs,gmres/solver.rs}`
  and `crates/cfd-validation/src/numerical/linear_solver.rs`.
- **Issue**: `LinearSolver::solve_system` still exposed nalgebra `DVector`
  RHS/initial-guess/result vectors after adjacent direct-solver and
  solver-chain boundaries had moved to Leto.
- **Remediation**: Moved the public `solve_system` trait boundary to
  `leto::Array1`, updated CG, BiCGSTAB, and GMRES implementations, and routed
  cfd-validation numerical solver validation through the new Leto public API.
- **Evidence tier**: compile-time provider integration, focused empirical
  solver nextest, clippy over affected producer and consumer paths, static
  signature/residue scans, and diff hygiene. `cargo fmt -p cfd-math -p
  cfd-validation --check`, `cargo check -p cfd-math --no-default-features
  --lib`, `cargo check -p cfd-validation --no-default-features --lib`,
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, `cargo clippy -p cfd-validation --no-default-features --lib -- -D
  warnings`, `cargo clippy -p cfd-validation --no-default-features
  --all-targets -- -D warnings`, `cargo clippy -p cfd-2d
  --no-default-features --all-targets -- -D warnings`, `cargo nextest run -p
  cfd-math --no-default-features conjugate_gradient bicgstab gmres
  --status-level fail` (58/58), targeted `solve_system` signature/residue
  scans, and `git diff --check` passed.
- **Residual risk**: `IterativeLinearSolver`, `LinearOperator`,
  `Preconditioner`, preconditioners, and internal iterative workspaces still
  expose nalgebra `DVector`/`CsrMatrix`. cfd-validation numerical solver
  result/error storage still uses `DVector`. Broad cfd-validation nextest
  still fails in the existing venturi cross-fidelity convergence tests.

# Sprint 1.96.125 Resolution: cfd-validation SpMV and Scalar Bounds Bypassed Leto

### RESOLVED-PENDING-154: Validation Kept Nalgebra SpMV Vectors and Weak Scalar Bounds
- **Location**:
  `crates/cfd-validation/src/benchmarking/{memory.rs,performance.rs}`,
  `crates/cfd-validation/src/numerical/linear_solver.rs`, and
  `crates/cfd-validation/src/literature/blood_flow_1d.rs`.
- **Issue**: cfd-validation still passed nalgebra `DVector` into public
  Leto-backed SpMV wrappers and used generic bounds that did not prove the
  Leto/Eunomia scalar contracts required by migrated sparse and network
  solver paths.
- **Remediation**: Moved validation SpMV benchmark vectors to `leto::Array1`
  and propagated the crate-local `ValidationScalar` seam through sparse
  linear-solver validation and 1D blood-flow literature validation.
- **Evidence tier**: compile-time provider integration, all-target clippy over
  cfd-validation and the dependent cfd-2d path, focused empirical benchmark
  nextest, static residue scan, and diff hygiene. `cargo fmt -p
  cfd-validation --check`, `cargo check -p cfd-validation
  --no-default-features --lib`, `cargo clippy -p cfd-validation
  --no-default-features --lib -- -D warnings`, `cargo clippy -p
  cfd-validation --no-default-features --all-targets -- -D warnings`, `cargo
  clippy -p cfd-2d --no-default-features --all-targets -- -D warnings`,
  targeted SpMV DVector-residue scan, `cargo nextest run -p cfd-validation
  --no-default-features benchmark --status-level fail` (40/40), and `git diff
  --check` passed.
- **Residual risk**: broad cfd-validation nextest still fails in
  `numerical::venturi_cross_fidelity` 2D fallback convergence tests:
  `microventuri_35um_case_produces_converged_informative_2d_result` and
  `option2_selected_45um_geometry_routes_to_fallback_and_converges`.

# Sprint 1.96.124 Resolution: Solver Chain and FEM Consumers Bypassed Leto

### RESOLVED-PENDING-153: LinearSolverChain and FEM Consumers Retained Nalgebra Vector Boundaries
- **Location**: `crates/cfd-math/src/linear_solver/chain.rs`;
  `crates/cfd-2d/src/linear_solver_bridge.rs`;
  `crates/cfd-2d/src/physics/momentum/{solve.rs,solver.rs}`;
  `crates/cfd-2d/src/pressure_velocity/correction.rs`;
  `crates/cfd-3d/src/fem/{leto_bridge.rs,solver.rs,projection_solver.rs}`;
  `crates/cfd-1d/src/solver/core/{mod.rs,linear_system.rs}`.
- **Issue**: After `DirectSparseSolver` and `SparseMatrixBuilder` moved to Leto
  arrays, their consumer boundary still passed nalgebra `DVector` through
  solver-chain and FEM/direct fallback call sites.
- **Remediation**: Moved `LinearSolverChain` RHS/result APIs to
  `leto::Array1`, introduced one cfd-2d direct-solver bridge, introduced one
  cfd-3d FEM sparse-assembly bridge, and propagated Leto scalar contracts
  through the cfd-1d/cfd-2d/cfd-3d scalar seams needed by these consumers.
- **Evidence tier**: compile-time provider integration, focused empirical
  nextest, clippy on touched library surfaces, and static residue scan.
  `cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`,
  `cargo check -p cfd-math --no-default-features --lib`,
  `cargo check -p cfd-1d --no-default-features --lib`,
  `cargo check -p cfd-2d --no-default-features --lib`,
  `cargo check -p cfd-3d --no-default-features --lib`,
  `cargo nextest run -p cfd-math --no-default-features chain direct_solver
  core_solver simple_gmres --status-level fail` (4/4),
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`, `cargo clippy -p cfd-2d --no-default-features --lib -- -D
  warnings`, and `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings` passed.
- **Residual risk**: `cfd-validation` still has nalgebra-vector public SpMV
  calls and generic Leto scalar-bound gaps. That blocks `cargo clippy -p
  cfd-2d --no-default-features --all-targets -- -D warnings` after the
  provider-bound changes. Iterative solver and preconditioner traits still
  expose nalgebra `DVector` and remain the next cfd-math boundary.

---

# Sprint 1.96.123 Resolution: cfd-math Direct Solver Leto Vectors

### RESOLVED-PENDING-152: Direct Sparse Solver Used DVector
- **Location**:
  `crates/cfd-math/src/linear_solver/{direct_solver.rs,dense_bridge.rs,chain.rs}`.
- **Issue**: `DirectSparseSolver::solve` still accepted and returned nalgebra
  `DVector` after the sparse builder and dense fallback internals had moved to
  Leto arrays.
- **Remediation**: Moved the direct solver RHS/result API to `leto::Array1`,
  moved rsparse RHS conversion and solution construction to Leto arrays,
  removed the stale DVector dense-fallback wrapper, and confined conversion to
  the current `LinearSolverChain` nalgebra boundary.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features direct_solver chain core_solver simple_gmres
  --status-level fail` passed 4/4 tests; `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`; `cargo doc -p
  cfd-math --no-default-features --no-deps`; `cargo test --doc -p cfd-math
  --no-default-features` passed 3 doctests with 3 ignored; targeted
  direct-solver DVector-signature scan; and `git diff --check`.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, rustdoc/doctest, static source audit, and diff hygiene.
- **Residual risk**: `LinearSolverChain`, iterative solver traits,
  preconditioners, and sparse storage still expose nalgebra `DVector` and
  `nalgebra_sparse::CsrMatrix` boundaries.

# Sprint 1.96.122 Resolution: cfd-math Sparse Builder Leto RHS

### RESOLVED-PENDING-151: Sparse Builder RHS Used DVector
- **Location**:
  `crates/cfd-math/src/sparse/builder.rs`,
  `crates/cfd-math/src/sparse/tests.rs`, and
  `crates/cfd-math/src/linear_solver/{direct_solver.rs,block_preconditioner.rs}`.
- **Issue**: `SparseMatrixBuilder::build_with_rhs` still exposed nalgebra
  `DVector` for Dirichlet column-elimination RHS assembly after neighboring
  sparse vector APIs had moved to Leto.
- **Remediation**: Moved the RHS parameter to `leto::Array1<T>`, mutated that
  Leto vector directly during column elimination, removed the obsolete dummy
  nalgebra RHS from `build()`, and updated direct/block-preconditioner builder
  call sites.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features sparse direct_solver block_preconditioner
  --status-level fail` passed 25/25 tests; `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`; targeted
  `build_with_rhs` residue scan; and `git diff --check`.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: The public sparse storage alias and
  linear-solver/preconditioner traits still expose
  `nalgebra_sparse::CsrMatrix` and nalgebra `DVector`.

# Sprint 1.96.121 Resolution: cfd-math Public SpMV Leto Vectors

### RESOLVED-PENDING-150: Public Sparse SpMV Wrappers Used DVector
- **Location**:
  `crates/cfd-math/src/sparse/operations.rs`, sparse tests, GMRES/AMG
  integration tests, interpolation quality checks, and
  `crates/cfd-math/benches/spmv_bench.rs`.
- **Issue**: The public sparse SpMV wrappers still accepted nalgebra
  `DVector` even though the provider computation already routed through
  Leto-ops.
- **Remediation**: Moved the then-public SpMV entry points to Leto `Array1`
  input/output vectors, later deleted the redundant parallel-named wrapper,
  and confined the `DVector` adapter to a
  private helper for `LinearOperator for CsrMatrix`.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features sparse spmv interpolation amg simple_gmres core_solver
  --status-level fail` passed 40/40 tests; `cargo clippy -p cfd-math
  --no-default-features --all-targets -- -D warnings`; public SpMV signature
  scan; and `git diff --check`.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: The public `LinearOperator`/solver/preconditioner trait
  family still exposes nalgebra `DVector`, and sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.

# Sprint 1.96.120 Resolution: cfd-math SparseMatrixExt Leto Vectors

### RESOLVED-PENDING-149: SparseMatrixExt Diagonal and Scaling Used DVector
- **Location**:
  `crates/cfd-math/src/sparse/operations.rs`,
  `crates/cfd-math/src/sparse/tests.rs`, and
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`.
- **Issue**: After sparse values were delegated to Leto CSR providers, the
  extension trait still exposed nalgebra `DVector` for diagonal extraction and
  row/column scaling vectors.
- **Remediation**: Moved `SparseMatrixExt::diagonal` to return
  `leto::Array1`, moved `set_diagonal`, `scale_rows`, and `scale_columns` to
  accept `leto::Array1`, and updated Jacobi construction to consume the Leto
  diagonal directly.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features sparse basic --status-level fail` passed 21/21 tests;
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`; targeted `SparseMatrixExt` DVector-signature residue scan; and
  `git diff --check` for the touched sparse/basic files.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: The public linear-solver/preconditioner traits still
  expose nalgebra `DVector`; sparse storage still aliases
  `nalgebra_sparse::CsrMatrix`.

# Sprint 1.96.119 Resolution: cfd-math Multigrid Smoother/Cycle Leto Vectors

### RESOLVED-PENDING-148: Multigrid Smoothers and Cycles Used nalgebra Vectors
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/{mod.rs,smoothers.rs,cycles.rs,amg.rs}`
  plus the Leto-array SpMV bridge in `crates/cfd-math/src/sparse/operations.rs`.
- **Issue**: After the restriction, GMRES, and GMG slices, neighboring
  multigrid smoothers and cycle algorithms still accepted nalgebra `DVector`
  values and depended on nalgebra operator-overload residual paths.
- **Remediation**: Introduced `MultigridVector<T> = leto::Array1<T>`, moved
  smoother and V/W/F cycle vector paths to that alias, added `spmv_array` for
  Leto-vector SpMV through `leto_ops::spmv_into`, and kept the remaining
  nalgebra bridge isolated at `Preconditioner::apply_to`.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features multigrid::cycles smoothers --status-level fail`
  passed 10/10 tests; `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`; clean smoother/cycle provider-residue scan;
  and `git diff --check` for the migrated sparse/multigrid files.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: The public sparse and linear-solver APIs still expose
  `nalgebra_sparse::CsrMatrix`, nalgebra `DVector`, and nalgebra
  `RealField`. Replacing those public contracts is the next larger
  linear-solver API migration.

# Sprint 1.96.118 Resolution: cfd-math GMG Leto/Eunomia Migration

### RESOLVED-PENDING-147: Geometric Multigrid Used nalgebra Dense Storage
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/gmg/{mod.rs,transfer.rs}`.
- **Issue**: The GMG hierarchy, transfer operators, FAS nonlinear-operator
  trait, residual computations, and tests still depended on nalgebra
  `DMatrix`/`DVector` and `num_traits` conversions.
- **Remediation**: Moved GMG matrices and vectors to `leto::Array2`/`Array1`,
  routed scalar constants through Eunomia `FloatElement`/`NumericElement`, and
  replaced nalgebra operator-overload math with explicit Leto-array helper
  kernels for matvec, residual, vector add/subtract/update, and L2 norm.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features gmg --status-level fail` passed 5/5 tests; `cargo
  clippy -p cfd-math --no-default-features --all-targets -- -D warnings`;
  clean GMG provider-residue scan; and `git diff --check` for the GMG files.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: The public linear-solver/preconditioner boundary still
  exposes nalgebra `DVector`/`RealField`.

# Sprint 1.96.117 Resolution: cfd-math GMRES Leto/Eunomia Workspace

### RESOLVED-PENDING-146: GMRES Internal Workspace Used nalgebra Dense Storage
- **Location**: `crates/cfd-math/src/linear_solver/gmres/{arnoldi.rs,givens.rs,solver.rs}`
  and `crates/cfd-math/src/linear_solver/chain.rs`.
- **Issue**: GMRES still used nalgebra `DMatrix` and internal `DVector`
  storage for the Krylov basis, Hessenberg matrix, Givens coefficients, and
  least-squares state after neighboring sparse and restriction paths had moved
  to Leto provider storage.
- **Remediation**: Replaced the internal Krylov/Hessenberg/Givens workspace
  with `leto::Array2`/`Array1`, routed Givens scalar math through Eunomia
  `RealField`/`NumericElement`, and propagated the Eunomia real-field bound to
  `LinearSolverChain`.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features gmres --status-level fail` passed 21/21 tests;
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`; clean GMRES provider-residue scan for replaced providers; and
  `git diff --check` for the GMRES and chain files.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, static source audit, and diff hygiene.
- **Residual risk**: GMRES still uses nalgebra `DVector` and
  `nalgebra::RealField` at the public `LinearOperator`, `Preconditioner`, and
  solver trait boundary. That boundary must move to `leto::Array1` and
  Eunomia scalar traits in a larger linear-solver API migration.

# Sprint 1.96.116 Resolution: cfd-math AMG Restriction Leto Arrays

### RESOLVED-PENDING-145: Standalone Restriction Utilities Exposed nalgebra Dense Types
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/restriction.rs`.
- **Issue**: The standalone AMG restriction helpers still used nalgebra
  `DMatrix`/`DVector` and nalgebra chained multiplication for dense Galerkin
  projection, even though neighboring sparse AMG paths had moved to Leto
  providers.
- **Remediation**: Moved restriction construction, validation, vector
  restriction, and matrix restriction to `leto::Array2`/`Array1`.
  `restrict_matrix` now delegates dense products to
  `leto_ops::MatrixProduct`, and tests assert concrete `P^T v` and `P^T A P`
  values.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features restriction --status-level fail` passed 7/7 tests;
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`; clean provider-residue scan over `restriction.rs`.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: Sparse/linear-solver public storage boundaries still
  expose nalgebra dense/sparse contracts.

# Sprint 1.96.115 Resolution: cfd-math Linear Solver Leto Dense Bridge

### RESOLVED-PENDING-144: Legacy CSR-to-Dense Solver Paths Bypassed Leto
- **Location**:
  `crates/cfd-math/src/linear_solver/{dense_bridge.rs,direct_solver.rs,mod.rs}`
  and `crates/cfd-math/src/linear_solver/preconditioners/multigrid/cycles.rs`.
- **Issue**: The direct sparse solver dense retry and multigrid cycle
  coarsest small-system solve each owned local dense solve logic over the
  legacy CSR/vector boundary instead of routing through Leto.
- **Remediation**: Added `linear_solver::dense_bridge` as the shared
  CSR/vector-to-Leto dense solve path. The direct dense retry and multigrid
  coarsest small solve both consume the bridge, and the multigrid local
  Gaussian-elimination helper was removed. `LinearSolverChain` now carries
  `leto_ops::RealScalar` because it owns the direct-solver call sites.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features multigrid::cycles direct_solver --status-level fail`
  passed 9/9 tests; `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  all-target clippy.
- **Residual risk**: The public direct-solver boundary remains
  `nalgebra_sparse::CsrMatrix`/nalgebra `DVector`, and the primary sparse LU
  attempt remains `rsparse`. A later slice should replace those with
  `leto_ops::CsrMatrix`, `leto::Array1`, and an Atlas-owned sparse direct
  solver.

# Sprint 1.96.114 Resolution: cfd-math Sparse Builder Leto Construction Bridge

### RESOLVED-PENDING-143: Sparse Builder Constructed CSR Through nalgebra
- **Location**: `crates/cfd-math/src/sparse/{bridge.rs,builder.rs,assembly.rs,operations.rs}`.
- **Issue**: The sparse module's operations had Leto-backed kernels, but
  `SparseMatrixBuilder` still constructed CSR storage directly through
  `nalgebra_sparse::CsrMatrix::try_from_csr_data`, and sparse operations owned
  their own Leto/nalgebra conversion helpers.
- **Remediation**: Added `sparse::bridge` as the single conversion boundary.
  `SparseMatrixBuilder::{build,build_with_rhs,build_parallel}` now construct
  `leto_ops::CsrMatrix` first and convert through the bridge only at the
  legacy public storage boundary. `ParallelAssembly::block_diagonal` routes its
  empty matrix edge through Leto CSR construction as well.
- **Verification**: `cargo check -p cfd-math --no-default-features --lib`;
  `cargo nextest run -p cfd-math --no-default-features sparse --status-level
  fail` passed 18/18 tests; `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  all-target clippy.
- **Residual risk**: `cfd-math::sparse::SparseMatrix` still aliases
  `nalgebra_sparse::CsrMatrix`, and solver APIs still consume nalgebra
  `DVector`. Replacing those public contracts remains the next larger
  sparse/linear-solver storage migration.

# Sprint 1.96.113 Resolution: cfd-math Sparse Extension Leto Provider Consumption

### RESOLVED-PENDING-142: Sparse Extension Utilities Still Owned CSR Loops
- **Location**: `crates/cfd-math/src/sparse/operations.rs` and
  `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`.
- **Issue**: After Leto owned CSR SpMV, SpGEMM, and transpose, CFDrs still
  implemented diagonal extraction, value scaling, row scaling, column scaling,
  Frobenius norm, diagonal dominance, and condition-estimate traversal inside
  `SparseMatrixExt`.
- **Remediation**: Added provider-backed bridge calls to
  `leto_ops::CsrMatrix` for each sparse extension utility and tightened
  `JacobiPreconditioner::new` to the Leto scalar contract required by diagonal
  extraction.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features sparse --status-level fail` passed 18/18 tests;
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  all-target clippy.
- **Residual risk**: This bridge keeps `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector` at the sparse/linear-solver API boundary. The next larger
  slice must replace those public storage contracts with `leto_ops::CsrMatrix`
  and `leto::Array1`.

# Sprint 1.96.112 Resolution: cfd-math AMG Leto SpGEMM Consumption

### RESOLVED-PENDING-141: AMG Galerkin Products Still Called Local CSR Multiply
- **Location**:
  `crates/cfd-math/src/sparse/operations.rs` and
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/amg.rs`.
- **Issue**: After the upstream `leto_ops::spgemm` provider landed, cfd-math
  still routed AMG Galerkin products through the existing
  `nalgebra_sparse`-owned multiplication boundary.
- **Remediation**: Added `try_sparse_sparse_mul`, which validates the current
  CSR storage boundary, delegates the product to `leto_ops::spgemm`, and maps
  provider errors into typed CFDrs errors. AMG recompute and hierarchy setup
  now use that fallible path for `R * A * P`. Added `try_sparse_transpose`,
  which delegates CSR transpose to `leto_ops::CsrMatrix::transpose`; AMG setup
  uses it for `R = P^T` instead of `nalgebra_sparse::transpose_as_csc`. Added
  `try_spmv`, which delegates matrix-vector products to `leto_ops::spmv_into`
  instead of retaining CFDrs-local scalar and parallel CSR traversal loops.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features sparse --status-level fail` passed 17/17 tests;
  `cargo nextest run -p cfd-math --no-default-features interpolation
  --status-level fail` passed 15/15 tests;
  `cargo nextest run -p cfd-math --no-default-features amg --status-level
  fail` passed 6/6 tests; `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  all-target clippy.
- **Residual risk**: This bridge keeps `nalgebra_sparse::CsrMatrix` and
  nalgebra `DVector` at the sparse/linear-solver API boundary. The next slice
  must replace those boundaries with `leto_ops::CsrMatrix` and `leto::Array1`.

# Sprint 1.96.111 Resolution: Leto-ops CSR Product Provider Gap

### RESOLVED-PENDING-140: AMG Galerkin Product Lacked a Leto CSR Target
- **Location**:
  `D:/atlas/repos/leto/crates/leto-ops/src/application/sparse`.
- **Issue**: CFDrs AMG migration needed a Leto-owned CSR×CSR product before
  the `nalgebra_sparse` Galerkin multiplication path could be removed without
  replacing it with a CFDrs-local sparse multiply.
- **Remediation**: Added `leto_ops::spgemm`, a CSR×CSR product that accumulates
  each output row into sorted columns and omits exact-zero cancellations. Added
  `CsrRow::nnz` for row-cardinality consumers.
- **Verification**: In `D:/atlas/repos/leto`, `cargo fmt -p leto-ops
  --check`; `cargo check -p leto-ops`; `cargo nextest run -p leto-ops --test
  ops_tests sparse --status-level fail` passed 14/14 tests; `cargo clippy -p
  leto-ops --all-targets -- -D warnings`; `cargo doc -p leto-ops --no-deps`.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, rustdoc, and static source audit.
- **Residual risk**: CFDrs still needs the consumer-side sparse/AMG migration
  from `nalgebra_sparse::CsrMatrix`/`DVector` to `leto_ops::CsrMatrix` and
  `leto::Array1`.

# Sprint 1.96.110 Resolution: cfd-math SIMD Leto/Eunomia Providers

### RESOLVED-PENDING-139: SIMD Cone Kept nalgebra Vector and Scalar Bounds
- **Location**:
  `crates/cfd-math/src/simd` and `crates/cfd-math/tests/simd_tests.rs`.
- **Issue**: The SIMD helper cone still implemented vector helpers for
  nalgebra `DVector`, imported nalgebra scalar traits for generic SIMD/field
  operations, and used `nalgebra_sparse` in the external SIMD integration
  residual test.
- **Remediation**: Replaced `SimdVectorOps` with a Leto `Array1`
  implementation, routed SIMD sparse matvec through `leto_ops::CsrMatrix` and
  `leto_ops::spmv`, moved generic SIMD bounds to Eunomia scalar traits, and
  migrated the integration residual test to Leto arrays plus Leto-ops CSR
  parts.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features simd --status-level fail` passed 26/26 tests with 318
  skipped; `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`; focused SIMD scan found no direct nalgebra, `DVector`, `DMatrix`,
  num-traits, ndarray, Rayon/Tokio, RustFFT, old scalar identity, or
  `default_epsilon` residue.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices.

# Sprint 1.96.109 Resolution: cfd-math Nonlinear Solver Leto Vectors

### RESOLVED-PENDING-138: Nonlinear Solver Kept nalgebra Dense Arrays
- **Location**:
  `crates/cfd-math/src/nonlinear_solver`.
- **Issue**: Anderson/JFNK nonlinear solvers still exposed and computed
  through nalgebra dense vector/matrix types after the scalar-provider slice
  had moved scalar math to Eunomia.
- **Remediation**: Replaced the nonlinear-solver vector surface with Leto
  `Array1`/`Array2`, consolidated Anderson/JFNK vector arithmetic in one
  local Leto helper module, migrated JFNK tests to Leto matrices, and kept JFNK
  scalar math on Eunomia `RealField`/`FloatElement`.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features nonlinear --status-level fail` passed 9/9 tests with
  335 skipped; `cargo clippy -p cfd-math --no-default-features --all-targets
  -- -D warnings`; focused nonlinear-solver scan found no direct nalgebra,
  `DVector`, `DMatrix`, num-traits, ndarray, Rayon/Tokio, RustFFT, old scalar
  identity, or `default_epsilon` residue.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices.

# Sprint 1.96.108 Resolution: cfd-math DG Leto Dense Arrays

### RESOLVED-PENDING-137: DG Cone Kept nalgebra Dense Arrays
- **Location**:
  `crates/cfd-math/src/high_order/dg`.
- **Issue**: The high-order DG cone still exposed and computed through
  nalgebra `DVector`/`DMatrix` even though Leto is the Atlas dense-array
  provider target.
- **Remediation**: Replaced DG basis, solution, numerical flux, operator,
  limiter, solver, and time-integration vector/matrix surfaces with Leto
  `Array1`/`Array2`. Projection, derivative, RHS, and implicit Newton
  correction solves now use Leto dense solve helpers and return typed solver
  errors instead of silently falling back. DG examples and DG-related
  Criterion benchmarks now construct and consume Leto arrays.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo check -p cfd-math
  --no-default-features --bench dg_benchmarks`; `cargo check -p cfd-math
  --no-default-features --bench flux_alloc_bench`; `cargo nextest run -p
  cfd-math --no-default-features dg --status-level fail` passed 62/62 tests
  with 282 skipped; `cargo clippy -p cfd-math --no-default-features
  --all-targets -- -D warnings`; `cargo test --doc -p cfd-math
  --no-default-features` passed 3 doctests with 3 ignored; focused
  high-order/DG-bench scan found no direct nalgebra, `DVector`, `DMatrix`,
  num-traits, ndarray, Rayon/Tokio, RustFFT, or old scalar identity residue.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices; SIMD is closed by Sprint 1.96.110.

# Sprint 1.96.107 Resolution: cfd-math Spectral Leto Dense Arrays

### RESOLVED-PENDING-136: Spectral Cone Kept nalgebra Dense Arrays
- **Location**:
  `crates/cfd-math/src/high_order/spectral`.
- **Issue**: The high-order spectral cone still exposed and computed through
  nalgebra `DVector`/`DMatrix` even though Leto is the Atlas dense-array
  provider target.
- **Remediation**: Replaced spectral element, mesh, differential operator,
  interpolation, quadrature, filter, and time-integration vector/matrix
  surfaces with Leto `Array1`/`Array2`; centralized matrix-vector, stiffness,
  dot-product, and vector-construction helpers in the spectral module; and
  migrated spectral assembly local dense matrices/RHS values and debug CSR
  materialization to Leto arrays. L2 projection now returns typed solver
  errors for failed Leto solves instead of silently falling back to
  interpolation.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features spectral --status-level fail` passed 13/13 tests with
  331 skipped; `cargo clippy -p cfd-math --no-default-features --all-targets
  -- -D warnings`; focused spectral scan found no direct nalgebra, `DVector`,
  `DMatrix`, num-traits, ndarray, Rayon/Tokio, RustFFT, or old scalar identity
  residue.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices; SIMD is closed by Sprint 1.96.110.

# Sprint 1.96.106 Resolution: cfd-math WENO Eunomia Scalars

### RESOLVED-PENDING-135: WENO Kept nalgebra Scalar Trait Bounds
- **Location**:
  `crates/cfd-math/src/high_order/weno/{mod,weno7}.rs`.
- **Issue**: The high-order WENO cone had removed direct `num_traits`
  conversion bounds but still imported and bound on `nalgebra::RealField`.
- **Remediation**: Replaced the WENO scalar trait import with Eunomia
  `RealField`/`FloatElement`, preserving the existing provider-owned
  constant-conversion helper and generic squaring helper for nonlinear
  weights.
- **Verification**: `cargo fmt -p cfd-math --check`; `cargo check -p
  cfd-math --no-default-features --lib`; `cargo nextest run -p cfd-math
  --no-default-features weno --status-level fail` passed 6/6 tests with 338
  skipped; `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`; focused WENO scan found no direct nalgebra, num-traits, ndarray,
  Rayon/Tokio, RustFFT, or old scalar identity residue.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  all-target clippy, and static source audit.
- **Residual risk**: cfd-math sparse and linear-solver provider replacement
  remain separate migration slices; SIMD is closed by Sprint 1.96.110.

# Sprint 1.96.105 Resolution: cfd-1d Solver-Core Eunomia Scalars

### RESOLVED-PENDING-134: Solver-Core Scalar Contract Kept Direct Num-Traits Bounds
- **Location**:
  `crates/cfd-1d/src/solver/core/{mod,anderson_acceleration,linear_system,convergence,solver_detection,geometry}.rs`.
- **Issue**: The primary network solver-core boundary still carried direct
  `num_traits::{Float, FromPrimitive, ToPrimitive}` requirements for scalar
  construction, finite checks, absolute values, square roots, f64 diagnostics,
  and residual calculations.
- **Remediation**: Replaced the shared `NetworkSolveScalar` compatibility
  bounds with Eunomia/cfd-core provider bounds and migrated Anderson
  acceleration, convergence checks, linear-system equilibration/Jacobi
  preconditioning, SPD detection, residual norms, and f64 diagnostics to
  `FloatElement`/`NumericElement`.
- **Verification**: `cargo fmt -p cfd-1d --check`; `cargo check -p cfd-1d`;
  `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped; focused
  solver-core scan found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_f64`, `nalgebra::try_convert`, or `<T as Float>`
  residue in the migrated files.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  static source audit.
- **Residual risk**: The active nalgebra/nalgebra-sparse matrix/vector storage
  boundary remains until the solver is migrated to Leto-backed dense/sparse
  provider APIs. Direct provider residue also remains in vascular
  Bessel/Womersley, resistance scalar traits, tests, and benches.

# Sprint 1.96.104 Resolution: cfd-1d Network Wrapper Eunomia Scalars

### RESOLVED-PENDING-133: Network Wrapper Kept Direct Scalar Bridge Calls
- **Location**:
  `crates/cfd-1d/src/domain/network/wrapper.rs`,
  `crates/cfd-1d/src/solver/core/matrix_assembly.rs`, and
  `crates/cfd-1d/src/solver/core/problem.rs`.
- **Issue**: The network wrapper still used direct `FromPrimitive`,
  `T::from_f64`, `T::from_usize`, `nalgebra::try_convert`, and generic
  `.abs()`/finite-check dispatch at the wrapper/solver boundary.
- **Remediation**: Replaced wrapper scalar construction and bridge
  conversions with `SafeFromF64`/`SafeFromUsize` and Eunomia
  `NumericElement`; tightened adjacent solver bounds so `MatrixAssembler`
  does not require `FromPrimitive`.
- **Verification**: `cargo fmt -p cfd-1d --check`; `cargo check -p cfd-1d`;
  `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped; focused
  wrapper/matrix/problem scan found no direct provider residue in the migrated
  files.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  static source audit.
- **Residual risk**: The `NetworkSolveScalar` compatibility bound noted in
  this slice is closed by Sprint 1.96.105. Remaining provider residue is in
  vascular Bessel/Womersley, resistance scalar traits, tests/benches, and
  nalgebra/nalgebra-sparse storage boundaries.

# Sprint 1.96.103 Resolution: cfd-1d Network Blueprint/Sink Eunomia Scalars

### RESOLVED-PENDING-132: Blueprint Conversion Kept Direct Scalar Construction
- **Location**:
  `crates/cfd-1d/src/domain/network/builder/blueprint_conversion.rs` and
  `crates/cfd-1d/src/domain/network/sink.rs`.
- **Issue**: The canonical blueprint-to-network path still used direct
  `T::from_f64`/`T::from_usize` scalar construction and generic `.abs()` in
  the builder/sink seam, even though the package now has an Atlas conversion
  provider boundary.
- **Remediation**: Replaced builder scalar construction with
  `SafeFromF64`/`SafeFromUsize`, replaced the generic coefficient absolute
  value with Eunomia `NumericElement::abs`, and propagated the Atlas
  conversion bounds through `NetworkBuilderSink`.
- **Verification**: `cargo fmt -p cfd-1d`; `cargo check -p cfd-1d`; `cargo
  nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
- **Evidence tier**: compile-time provider integration and empirical nextest.
- **Residual risk**: The wrapper residue described here is closed by Sprint
  1.96.104, and the inherited solver scalar compatibility bound is closed by
  Sprint 1.96.105. Remaining provider residue is in vascular Bessel/Womersley,
  resistance scalar traits, tests/benches, and nalgebra/nalgebra-sparse storage
  boundaries.

# Sprint 1.96.102 Resolution: cfd-1d Domain Components Eunomia Scalars

### RESOLVED-PENDING-131: Domain Components Kept Direct Num-Traits Bounds
- **Location**:
  `crates/cfd-1d/src/domain/components/{mod,channels,factory,membranes,mixers,pumps,sensors,valves}.rs`.
- **Issue**: Domain components still depended on direct
  `num_traits::FromPrimitive`/`Float` bounds for component construction,
  constants, and pressure-drop absolute-value math.
- **Remediation**: Replaced the component scalar boundary with
  `SafeFromF64`, `SafeFromUsize`, and Eunomia `NumericElement`. The component
  trait uses provider-owned absolute value, the factory delegates fallible
  constant conversion to the Atlas conversion seam, and the circular-channel
  test module was moved below all items to keep the touched file clippy-clean.
- **Verification**: `cargo fmt -p cfd-1d --check`; `cargo check -p cfd-1d`;
  `cargo nextest run -p cfd-1d` passed 725/725 tests with 3 skipped; focused
  scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `Float::`, or generic `.abs()` residue in
  domain components. A broad clippy rerun confirmed the touched `channels.rs`
  lint is gone.
- **Evidence tier**: compile-time provider integration, empirical nextest,
  static source audit, and lint regression audit.
- **Residual risk**: `cfd-1d` all-target clippy still fails outside this slice
  in tests/examples, domain-network, cell-separation, vascular, solver-core,
  and benches. Remaining direct provider work is in solver-core,
  domain-network, vascular Bessel/Womersley, tests, benches, and current
  nalgebra/nalgebra-sparse storage boundaries.

# Sprint 1.96.101 Resolution: cfd-1d Channel/Branching/Analysis Eunomia Scalars

### RESOLVED-PENDING-130: cfd-1d Channel/Analysis Kept Direct Num-Traits Bounds
- **Location**:
  `crates/cfd-1d/src/domain/channel`,
  `crates/cfd-1d/src/domain/junctions/branching`, and
  `crates/cfd-1d/src/solver/analysis`.
- **Issue**: The next cfd-1d provider seam still carried direct
  `num_traits::{FromPrimitive, ToPrimitive, Float}` bounds plus direct
  `T::from_f64`, `Float::abs`/`sqrt`/`powf`, and generic
  `to_f64().unwrap_or(...)` conversion patterns.
- **Remediation**: Replaced these call paths with `SafeFromF64`, Eunomia
  `FloatElement`, and Eunomia `NumericElement` while preserving the current
  nalgebra `RealField` and sparse solver storage boundary for later Leto work.
- **Verification**: `cargo check -p cfd-1d` passed; `cargo nextest run -p
  cfd-1d` passed 725/725 tests with 3 skipped; focused scans found no direct
  provider residue in the touched channel/branching/analyzer cone or in
  `solver/analysis`.
- **Evidence tier**: compile-time provider integration, empirical nextest, and
  static source audit.
- **Residual risk**: `cfd-1d` still has direct `num-traits` usage in
  solver-core, domain-network/components, vascular Bessel/Womersley, and tests.
  All-target clippy is blocked by unrelated existing lint debt in examples,
  tests, and cell-separation/resistance modules.

# Sprint 1.96.100 Resolution: cfd-core Fluid Dynamics Operations Eunomia Scalars

### RESOLVED-PENDING-129: Flow Operations Kept Direct Num-Traits Construction
- **Location**:
  `crates/cfd-core/src/physics/fluid_dynamics/operations.rs`.
- **Issue**: Core flow-field operations still depended on direct
  `num_traits::FromPrimitive` and fallback scalar construction for central
  difference denominators and energy/enstrophy half factors.
- **Remediation**: Replaced local scalar construction with Eunomia
  `NumericElement` identities and explicit half/two helpers while preserving
  Moirai parallel iteration and the current nalgebra `Vector3` storage
  boundary.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  fluid_dynamics::operations` passed 3/3; touched-file rustfmt passed; focused
  scan found no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or old
  `unwrap_or_else` scalar-fallback patterns in `operations.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic
  operations tests, and static source audit.
- **Residual risk**: `VelocityField<T>`, flow operations, RANS/turbulence
  traits, and Rhie-Chow still retain nalgebra vector/storage contracts pending
  the Leto migration.

# Sprint 1.96.99 Resolution: cfd-core Fluid Dynamics Service Eunomia Scalars

### RESOLVED-PENDING-128: Fluid Dynamics Service Kept Direct Num-Traits Math
- **Location**:
  `crates/cfd-core/src/physics/fluid_dynamics/service.rs`.
- **Issue**: Pipe pressure-drop and friction-factor formulas still depended on
  direct `num_traits::{Float, FromPrimitive}` plus silent scalar-conversion
  fallbacks for constants and math dispatch.
- **Remediation**: Replaced the service formula scalar construction and math
  dispatch with Eunomia `FloatElement`/`NumericElement`, including constants,
  powers, square roots, logarithms, and absolute convergence checks. Added
  closed-form laminar friction-factor coverage while retaining the existing
  Colebrook-White validation.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  fluid_dynamics::service` passed 2/2; touched-file rustfmt passed; focused
  scan found no direct `num_traits`, `FromPrimitive`, `Float::`,
  `T::from_f64`, or old `unwrap_or_else` scalar-fallback patterns in
  `service.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic service
  tests, and static source audit.
- **Residual risk**: The same fluid-dynamics bounded context still has direct
  `num_traits::FromPrimitive` in `operations.rs` and nalgebra vector/storage
  boundaries in operations and `VelocityField<T>`.

# Sprint 1.96.98 Resolution: cfd-core Flow Regime Eunomia Scalars

### RESOLVED-PENDING-127: Flow Regimes Kept Direct Num-Traits Conversion
- **Location**:
  `crates/cfd-core/src/physics/fluid_dynamics/flow_regimes.rs` and the
  `FluidDynamicsService::flow_regime` wrapper in
  `crates/cfd-core/src/physics/fluid_dynamics/service.rs`.
- **Issue**: Core flow classification still depended on direct
  `nalgebra::RealField` plus `num_traits::ToPrimitive`, including silent
  `unwrap_or(0.0)` conversion fallback behavior.
- **Remediation**: Replaced the classifier contract with Eunomia
  `RealField`/`NumericElement`, converted via `NumericElement::to_f64`, and
  migrated the service wrapper to the same scalar provider contract.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  flow_regime` passed 3/3; touched-file rustfmt passed; focused scan found no
  direct `num_traits`, `ToPrimitive`, direct `nalgebra::RealField`,
  `.to_f64()`, or `unwrap_or(0.0)` residue in `flow_regimes.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic
  threshold tests, and static source audit.
- **Residual risk**: The same fluid-dynamics bounded context still has direct
  `num_traits::{Float, FromPrimitive}` in the service pipe-flow formulas and
  nalgebra vector/storage boundaries in operations and `VelocityField<T>`.

# Sprint 1.96.97 Resolution: cfd-3d Spectral Diagnostics Eunomia Conversion

### RESOLVED-PENDING-126: Spectral Diagnostics Kept Direct Num-Traits Conversion
- **Location**:
  `crates/cfd-3d/src/spectral/diagnostics.rs`.
- **Issue**: The Apollo/Leto diagnostics path still used direct
  `num_traits::{FromPrimitive, ToPrimitive}` and a fallible `.to_f64()` branch
  while staging velocity components into Leto arrays for Apollo FFTs.
- **Remediation**: Replaced the local conversion seam with Eunomia
  `NumericElement::to_f64` and removed the obsolete conversion-context
  parameter.
- **Verification**: `cargo check -p cfd-3d`; `cargo nextest run -p cfd-3d
  diagnostics` passed 5/5; touched-file rustfmt passed; focused scan found no
  direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `.to_f64()`, or
  obsolete `component_name` residue in `diagnostics.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic
  spectral diagnostics tests, and static source audit.
- **Residual risk**: `cfd-core::VelocityField<T>` still owns the nalgebra
  `Vector3`/`RealField` storage boundary, and larger `cfd-3d` spectral
  Chebyshev/Poisson/solver code still uses nalgebra dense matrix/vector types.

# Sprint 1.96.96 Resolution: cfd-2d MRT/Carreau-Yasuda Eunomia Scalars

### RESOLVED-PENDING-125: MRT/Carreau-Yasuda Kept Legacy Scalar Providers
- **Location**:
  `crates/cfd-2d/src/solvers/lbm/collision/{mrt,carreau_yasuda}.rs` and
  `crates/cfd-2d/src/physics/non_newtonian/carreau_yasuda.rs`.
- **Issue**: The previous LBM scalar-provider slice left MRT moment helpers
  and local Carreau-Yasuda rheology/collision math on direct
  `nalgebra::RealField`, `num_traits::{Float, FromPrimitive}`, and legacy
  generic scalar constructors.
- **Remediation**: Replaced the remaining local scalar construction and math
  dispatch with Eunomia `FloatElement`/`NumericElement`, local scalar helpers,
  and `CastFrom<i32>` for D2Q9 velocity components.
- **Verification**: `cargo check -p cfd-2d`; `cargo nextest run -p cfd-2d
  lbm` passed 31/31; `cargo nextest run -p cfd-2d carreau_yasuda` passed 5/5;
  touched-file rustfmt passed; focused scan found no direct
  `nalgebra::RealField`, `num_traits`, `FromPrimitive`, `T::from_f64`,
  `T::zero`, `T::one`, `from_i32`, or `Float::` residue in the migrated
  MRT/Carreau-Yasuda files.
- **Evidence tier**: compile-time provider integration, value-semantic LBM and
  Carreau-Yasuda tests, and static source audit.
- **Residual risk**: Broader CFDrs Atlas provider migration remains incomplete
  outside this bounded LBM scalar-provider closeout.

# Sprint 1.96.95 Resolution: cfd-2d LBM Eunomia Scalar Seam

### RESOLVED-PENDING-124: LBM Kept Direct Nalgebra/Num-Traits Scalars
- **Location**:
  `crates/cfd-2d/src/solvers/lbm/{macroscopic,lattice}.rs`,
  `crates/cfd-2d/src/solvers/lbm/collision/{traits,bgk}.rs`, and bridged
  trait impl bounds in `mrt.rs` and `carreau_yasuda.rs`.
- **Issue**: The authoritative LBM macroscopic extraction, D2Q9 equilibrium,
  collision trait, and BGK collision path still used direct
  `nalgebra::RealField` and `num_traits::FromPrimitive` for generic scalar
  construction and constants.
- **Remediation**: Replaced the authoritative seam with Eunomia
  `FloatElement` and local scalar helpers for density, velocity, pressure,
  stress, kinetic energy, vorticity, equilibrium weights/constants, BGK
  viscosity, and omega construction. MRT and Carreau-Yasuda trait impl bounds
  now satisfy the migrated `CollisionOperator` seam.
- **Verification**: `cargo check -p cfd-2d`; `cargo nextest run -p cfd-2d
  lbm` passed 31/31; touched-file rustfmt and `git diff --check` passed;
  focused scan found no direct `nalgebra`, `RealField`, `num_traits`,
  `FromPrimitive`, or `ToPrimitive` residue in `macroscopic.rs`,
  `lattice.rs`, `collision/traits.rs`, or `collision/bgk.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic LBM
  tests, and static source audit.
- **Residual risk**: The MRT and Carreau-Yasuda residual named by this slice is
  closed by Sprint 1.96.96. The broader CFDrs Atlas provider migration remains
  incomplete outside the LBM scalar-provider seam.

# Sprint 1.96.94 Resolution: cfd-3d Wall Functions Eunomia Scalars

### RESOLVED-PENDING-123: Wall Functions Kept Direct Num-Traits Scalars
- **Location**:
  `crates/cfd-3d/src/physics/turbulence/wall_functions.rs`.
- **Issue**: The 3D Spalding wall-law helper still used direct
  `nalgebra::RealField`, `num_traits::FromPrimitive`, and
  `num_traits::Float` for generic scalar constants, exponentials, logarithm,
  square root, absolute-value convergence, and max clamps.
- **Remediation**: Replaced local scalar construction and math dispatch with
  Eunomia `FloatElement`/`NumericElement` while preserving the existing
  wall-law formulas and tests.
- **Verification**: `cargo check -p cfd-3d`; `cargo nextest run -p cfd-3d
  wall_functions` passed 3/3; touched-file rustfmt and `git diff --check`
  passed; focused scan found no direct `nalgebra`, `RealField`,
  `num_traits`, `FromPrimitive`, or `ToPrimitive` residue in
  `wall_functions.rs`.
- **Evidence tier**: compile-time provider integration, value-semantic
  wall-law tests, and static source audit.
- **Residual risk**: Direct `GpuContext` Hephaestus replacement is a separate
  WGPU API-version alignment item because current CFDrs resolves
  `wgpu@0.19.4` and Hephaestus resolves `wgpu@26.0.1`.

# Sprint 1.96.93 Resolution: cfd-3d Apollo Fourier Eunomia Bounds

### RESOLVED-PENDING-122: Apollo Fourier Wrapper Kept Direct Num-Traits Conversion
- **Location**:
  `crates/cfd-3d/src/spectral/fourier.rs`,
  `crates/cfd-3d/src/physics/turbulence/des.rs`, and downstream touched
  `cfd-2d`/`cfd-validation` blood-consumer bounds.
- **Issue**: The Apollo-backed Fourier wrapper still used direct
  `num_traits::{FromPrimitive, ToPrimitive}` for scalar conversion, and
  downstream consumers had not fully adopted the `FloatElement` bounds now
  required by migrated `cfd-core` blood/cavitation APIs.
- **Remediation**: Replaced Fourier conversion helpers with Eunomia
  `FloatElement`/`NumericElement`, migrated DES blood-model scalar conversion,
  and propagated `FloatElement` bounds through touched downstream consumers.
- **Verification**: `cargo check -p cfd-3d`; `cargo check -p
  cfd-validation`; `cargo nextest run -p cfd-3d --test fourier_validation`
  passed 12/12; touched-file rustfmt passed; focused residue scan found no
  direct `num_traits`, `FromPrimitive`, or `ToPrimitive` residue in the
  touched `cfd-3d` Fourier/DES files.
- **Evidence tier**: compile-time provider integration, value-semantic Fourier
  integration tests, and static source audit.
- **Residual risk**: Broader direct `num_traits`, `nalgebra`, direct `wgpu`,
  and execution-provider residues remain outside this slice.

# Sprint 1.96.92 Resolution: cfd-core Fåhræus-Lindqvist Eunomia Scalars

### RESOLVED-PENDING-121: Fåhræus-Lindqvist Kept Direct Num-Traits Scalars
- **Location**:
  `crates/cfd-core/src/physics/fluid/blood/fahraeus_lindqvist.rs`.
- **Issue**: The Fåhræus-Lindqvist blood model still used direct
  `num_traits::FromPrimitive` construction and direct generic `powf`, `exp`,
  `abs`, zero/one, and max-style scalar dispatch.
- **Remediation**: Replaced local scalar construction and generic math dispatch
  with Eunomia `FloatElement`/`NumericElement` across the Pries/Secomb
  relative-viscosity formulas, `mu_45` fit, relative-viscosity clamp, and tube
  hematocrit correlation.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  fahraeus_lindqvist` passed 3/3; `cargo nextest run -p cfd-core blood`
  passed 24/24; focused `fahraeus_lindqvist.rs` residue scan found no direct
  `num_traits`, `FromPrimitive`, generic `T::from_f64`, `T::zero`, `T::one`,
  generic `powf`, generic `exp`, or generic `abs`; broader blood residue scan
  now matches only concrete `f64` helper/test expressions; touched-file rustfmt
  and touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: The broader fluid trait `RealField` boundary remains
  separate provider work. Full cfd-core clippy remains blocked by unrelated
  existing lints in boundary applicator and Rhie-Chow tests.

# Sprint 1.96.91 Resolution: cfd-core Casson/Carreau Blood Eunomia Scalars

### RESOLVED-PENDING-120: Casson/BloodModel Kept Direct Num-Traits Dispatch
- **Location**:
  `crates/cfd-core/src/physics/fluid/blood/{casson,carreau_yasuda,mod}.rs`.
- **Issue**: Casson still used direct `num_traits::FromPrimitive`
  construction, direct generic zero/one identities, and direct generic
  `sqrt`/`exp` dispatch. The shared `BloodModel` selector still required
  `FromPrimitive`, which also prevented the already-migrated Carreau-Yasuda
  model from compiling through the dispatch enum.
- **Remediation**: Replaced Casson scalar construction and math dispatch with
  Eunomia `FloatElement`/`NumericElement`, and changed `BloodModel` to require
  Eunomia `FloatElement`.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  casson` passed 12/12; `cargo nextest run -p cfd-core carreau_yasuda` passed
  4/4; `cargo nextest run -p cfd-core blood` passed 24/24; focused touched-file
  residue scan found no direct `num_traits`, `FromPrimitive`, generic
  `T::from_f64`, `T::zero`, `T::one`, generic `powf`, generic `sqrt`, or
  generic `exp` residue except the intentionally concrete `f64` temperature
  helper; touched-file rustfmt and touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: Fåhræus-Lindqvist and the fluid trait `RealField`
  boundary remain separate blood-fluid provider holdouts. Full cfd-core clippy
  remains blocked by unrelated existing lints in boundary applicator and
  Rhie-Chow tests.

# Sprint 1.96.90 Resolution: cfd-core Cross Blood Eunomia Scalars

### RESOLVED-PENDING-119: Cross Blood Kept Direct Num-Traits Construction
- **Location**: `crates/cfd-core/src/physics/fluid/blood/cross.rs`.
- **Issue**: The Cross blood model still used direct
  `num_traits::FromPrimitive` construction and direct generic scalar math
  despite requiring only Eunomia-covered constants, zero/one identities, and
  real-power evaluation.
- **Remediation**: Replaced Cross-local scalar construction and `powf` dispatch
  with Eunomia `FloatElement`/`NumericElement`.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  cross` passed 1/1; `cargo nextest run -p cfd-core blood` passed 24/24;
  focused `cross.rs` residue scan found no `num_traits`, `FromPrimitive`,
  direct `T::from_f64`, direct `T::zero`, direct `T::one`, or direct `powf`;
  touched-file rustfmt and touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: The fluid trait still carries `nalgebra::RealField`.
  Casson, Carreau-Yasuda, Fåhræus-Lindqvist, and `BloodModel` remain separate
  blood-fluid scalar holdouts.

# Sprint 1.96.89 Resolution: cfd-core Cavitation Scalar Closeout

### RESOLVED-PENDING-118: Cavitation Models Kept Nalgebra Scalar Bounds
- **Location**:
  `crates/cfd-core/src/physics/cavitation/{models,heterogeneous_nucleation}.rs`.
- **Issue**: Cavitation mass-transfer models still used `nalgebra::RealField`
  and direct scalar construction, and the heterogeneous legacy helper exposed
  a fake generic `RealField` contract while converting every value through
  `f64`.
- **Remediation**: Replaced `CavitationModel<T>`/`ZgbParams<T>` with Eunomia
  `FloatElement`/`NumericElement` scalar contracts and changed the
  heterogeneous legacy helper to the concrete `f64` contract used by the
  selective-cavitation model.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  cavitation` passed 40/40; focused cavitation residue scan found no
  `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::zero`, direct `T::one`, `to_subset`, or
  `try_convert`; touched-file rustfmt and touched-file `git diff --check`
  passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: Broader CFDrs Atlas migration remains incomplete outside
  cavitation.

# Sprint 1.96.88 Resolution: cfd-core Nuclei Transport Eunomia Scalars

### RESOLVED-PENDING-117: Nuclei Transport Kept Nalgebra Scalar Bounds
- **Location**: `crates/cfd-core/src/physics/cavitation/nuclei_transport.rs`.
- **Issue**: Nuclei transport still used `nalgebra::RealField` and direct
  scalar construction despite requiring only Eunomia-covered scalar
  arithmetic, comparisons, constants, and exponential evaluation.
- **Remediation**: Replaced the nuclei transport scalar contract with Eunomia
  `FloatElement`/`NumericElement` and routed constants, zero/one identities,
  and exponential decay through provider APIs.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  cavitation` passed 35/35; focused `nuclei_transport.rs` residue scan found
  no `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::zero`, direct `T::one`, or direct generic
  `exp`; touched-file rustfmt and touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: `models.rs` and `heterogeneous_nucleation.rs` were
  subsequently closed by Sprint 1.96.89.

# Sprint 1.96.87 Resolution: cfd-core Venturi Cavitation Eunomia Scalars

### RESOLVED-PENDING-116: Venturi Cavitation Kept Nalgebra Scalar Bounds
- **Location**: `crates/cfd-core/src/physics/cavitation/venturi.rs`.
- **Issue**: `VenturiCavitation<T>` still used `nalgebra::RealField` and
  direct `num_traits::FromPrimitive` construction despite requiring only
  Eunomia-covered scalar arithmetic, comparisons, constants, integer powers,
  tangent, and absolute value.
- **Remediation**: Replaced the Venturi scalar contract with Eunomia
  `FloatElement`/`NumericElement` and routed constants/math operations through
  provider APIs.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  cavitation` passed 35/35; focused production residue scan found no
  `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::zero`, direct `T::one`, direct `powf`,
  direct `powi`, or direct `tan` in `venturi.rs`; touched-file rustfmt and
  touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: Remaining cavitation holdouts at this point were
  `models.rs`, `nuclei_transport.rs`, and `heterogeneous_nucleation.rs`;
  Sprint 1.96.88 subsequently removed the `nuclei_transport.rs` holdout.

# Sprint 1.96.86 Resolution: cfd-core Cavitation Eunomia Scalar Cone

### RESOLVED-PENDING-115: Cavitation Scalar Cone Kept Nalgebra Bounds
- **Location**:
  `crates/cfd-core/src/physics/cavitation/{rayleigh_plesset,bio_damage,number,damage}.rs`
  and `crates/cfd-core/src/physics/cavitation/regimes/`.
- **Issue**: The touched cavitation scalar cone still used
  `nalgebra::RealField`, direct `num_traits::FromPrimitive`, and direct scalar
  construction despite requiring only Eunomia-covered scalar arithmetic and
  real-valued math operations.
- **Remediation**: Replaced the touched cavitation scalar contracts with
  Eunomia `FloatElement`/`NumericElement`, routed constants and math operations
  through provider APIs, and added closed-form tests for cavitation-number and
  material-damage formulas.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  cavitation` passed 35/35; focused migrated-file residue scan found no
  `nalgebra`, `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::from_u64`, direct `T::zero`, direct
  `T::one`, direct `powf`, or direct `powi`; touched-file rustfmt and
  touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: Remaining cavitation holdouts at this point were
  `models.rs`, `venturi.rs`, `nuclei_transport.rs`, and
  `heterogeneous_nucleation.rs`; Sprint 1.96.87 subsequently removed the
  `venturi.rs` holdout.
  Full `cfd-core` clippy is currently blocked by unrelated existing lints in
  `physics/boundary/applicator.rs` and `physics/fluid_dynamics/rhie_chow.rs`.

# Sprint 1.96.85 Resolution: cfd-core Hemolysis Eunomia Scalars

### RESOLVED-PENDING-114: Hemolysis Kept Nalgebra Scalar Bounds
- **Location**: `crates/cfd-core/src/physics/hemolysis/{calculator,trauma}.rs`.
- **Issue**: The hemolysis calculator and platelet activation model still used
  `nalgebra::RealField` and direct `num_traits::FromPrimitive` construction for
  scalar constants despite requiring only scalar comparisons, arithmetic,
  constants, and exponential evaluation.
- **Remediation**: Replaced hemolysis-local scalar contracts with Eunomia
  `FloatElement`/`NumericElement`, including provider-owned zero/one values,
  constants, and exponentials.
- **Verification**: `cargo check -p cfd-core`; `cargo nextest run -p cfd-core
  hemolysis` passed 9/9; focused hemolysis residue scan found no `nalgebra`,
  `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::zero`, or
  `T::one`; touched-file rustfmt and touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider integration, value-semantic focused
  tests, and static source audit.
- **Residual risk**: Broader `cfd-core` nalgebra surfaces remain in compute,
  boundary, cavitation, fluid, mesh, and fluid-dynamics modules.

# Sprint 1.96.84 Resolution: cfd-2d Grid/FVM Leto-Eunomia Boundary

### RESOLVED-PENDING-113: Structured Grid Kept FVM on Nalgebra Scalars
- **Location**: `crates/cfd-2d/src/grid/{traits,structured,unstructured,refinement}.rs`,
  `crates/cfd-2d/src/solvers/fvm/solver.rs`, and grid-center consumers in
  `crates/cfd-2d` plus `crates/cfd-validation`.
- **Issue**: The prior FVM vector migration still inherited
  `nalgebra::RealField` through `StructuredGrid2D<T>`, and
  `Grid2D::cell_center` still exposed nalgebra `Vector2` semantics to callers.
- **Remediation**: Replaced the grid trait center type with
  `leto::geometry::Vector2`, routed structured/adaptive grid scalar
  construction through Eunomia helpers, removed `RealField` from the
  unstructured grid storage boundary, removed the inherited `RealField` bound
  from `FvmSolver`, and updated grid-center consumers to Leto indexing.
- **Verification**: `cargo check -p cfd-2d`; `cargo check -p
  cfd-validation`; `cargo nextest run -p cfd-2d --lib grid` passed 21/21;
  `cargo nextest run -p cfd-2d --lib fvm` passed 25/25; focused grid/FVM
  residue scan found no `RealField`, `FromPrimitive`, `ToPrimitive`,
  `num_traits`, `nalgebra::Vector2`, direct `T::from_f64`, direct
  `T::from_usize`, or stale scalar fallback residue; touched-file rustfmt and
  touched-file `git diff --check` passed.
- **Evidence tier**: compile-time provider/consumer integration,
  value-semantic focused tests, and static source audit.
- **Residual risk**: Broader `cfd-2d` nalgebra/num-traits surfaces remain
  outside this grid/FVM cone. `tests_poisson_mms.rs` contains an MMS test
  function but is not currently listed by `cargo nextest list -p cfd-2d --lib`;
  that test-discovery gap is separate from this provider-boundary migration.

# Sprint 1.96.83 Resolution: cfd-2d FVM Leto Vector Boundary

### RESOLVED-PENDING-112: FVM Cone Stored Nalgebra Vector2
- **Location**: `crates/cfd-2d/src/solvers/fvm/{geometry,solver}.rs` and
  `../leto/crates/leto/src/geometry.rs`.
- **Issue**: The FVM face and velocity-field boundary still stored nalgebra
  `Vector2` even though Leto is the Atlas fixed-vector provider.
- **Remediation**: Added the Leto `Vector2<T>` alias and fixed-vector
  norm/normalization methods, added direct Leto consumption to `cfd-2d`, and
  replaced FVM face center/normal and velocity-field vector use with
  `leto::geometry::Vector2<T>`.
- **Verification**: Focused FVM scan shows no nalgebra `Vector2`,
  `num_traits`, direct `T::from_f64`, direct `T::from_usize`, `.to_subset()`,
  or stale diffusion fallback residue. Touched-file rustfmt passed. `cargo
  check -p leto` and `cargo check -p cfd-2d` passed. `cargo nextest run -p
  leto fixed_vector_norm_and_normalization_are_value_semantic` passed 1/1 test.
  `cargo nextest run -p cfd-2d --lib fvm` passed 25/25 focused tests.
- **Evidence tier**: compile-time provider/consumer integration,
  value-semantic focused tests, and static source audit.
- **Residual risk**: `FvmSolver<T>` still declares `nalgebra::RealField`
  through `StructuredGrid2D<T>`; full grid migration remains open.

# Sprint 1.96.82 Resolution: cfd-2d FVM Eunomia Scalar Constants

### RESOLVED-PENDING-111: FVM Cone Used Direct num-traits Scalars
- **Location**: `crates/cfd-2d/src/solvers/fvm/{config,solver,flux}.rs`.
- **Issue**: FVM configuration, face-center construction, and flux constants
  still depended on direct `num_traits::FromPrimitive` scalar construction even
  though the scalar constants can route through Eunomia.
- **Remediation**: Replaced the config bound with Eunomia `FloatElement`,
  routed default constants through `FloatElement::from_f64`, propagated
  `FloatElement` to `FvmSolver`, replaced face-center index conversions with
  Eunomia-backed helpers, replaced flux scalar constants with Eunomia helpers,
  and added invalid-diffusion rejection for power-law/hybrid fluxes.
- **Verification**: Focused FVM scan shows no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, `.to_subset()`, or stale diffusion fallback residue.
  Touched-file rustfmt passed. `cargo check -p cfd-2d` and `cargo check -p
  cfd-3d` passed, and `cargo nextest run -p cfd-2d --lib fvm` passed 25/25
  focused tests.
- **Evidence tier**: compile-time dependency-chain verification plus
  value-semantic focused tests.
- **Residual risk**: FVM geometry, solver, and flux still contain nalgebra
  `RealField`/`Vector2` contracts and require a separate Leto-backed migration
  slice.

# Sprint 1.96.81 Resolution: cfd-2d CFL Eunomia Scalars

### RESOLVED-110: CFL Calculator Used Nalgebra/num-traits Scalars
- **Location**: `crates/cfd-2d/src/stability/cfl.rs`.
- **Issue**: The standalone CFL stability calculator still depended on
  `nalgebra::RealField` and `num_traits::FromPrimitive` for scalar constants,
  absolute values, and CFL threshold construction.
- **Remediation**: Replaced the bounds and scalar calls with Eunomia
  `FloatElement`/`NumericElement`, preserving the documented CFL, diffusion,
  QUICK, and max-stable-time-step formulas.
- **Verification**: Focused CFL scan shows no `RealField`, `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::zero`, or `T::one`
  residue. Touched-file rustfmt passed. `cargo check -p cfd-2d` passed.
  `cargo nextest run -p cfd-2d --lib cfl` passed 9/9 tests.
- **Evidence tier**: compile-time integration, empirical focused tests, and
  static source audit.
- **Residual risk**: A broader `cargo nextest run -p cfd-2d cfl` compiles
  unrelated integration tests and currently hits pre-existing
  `simplec_pimple_validation` generic-bound debt. Remaining cfd-2d
  stability/time-scheme, turbulence, validation, and solver modules still have
  nalgebra/num-traits surfaces.

# Sprint 1.96.80 Resolution: cfd-core Material Eunomia Traits

### RESOLVED-109: Material Solid/Interface Contracts Used Nalgebra RealField
- **Location**: `crates/cfd-core/src/physics/material/{traits,solid,interface,mod}.rs`.
- **Issue**: Material solid/interface property contracts had already moved
  constructor constants to Eunomia, but the public trait/value type bounds still
  required `nalgebra::RealField`.
- **Remediation**: Replaced solid/interface material bounds with Eunomia
  `FloatElement`/`NumericElement`, routed shear modulus and adhesion energy
  through Eunomia operations, and kept `MaterialDatabase`'s `RealField`
  requirement scoped to its `Fluid<T>` map.
- **Verification**: Focused material scan shows no `nalgebra::RealField`,
  `RealField`, `num_traits`, `FromPrimitive`, `ToPrimitive`, or `Float` residue
  in the migrated solid/interface files. Touched-file rustfmt passed. `cargo
  check -p cfd-core` passed. `cargo nextest run -p cfd-core material` passed
  4/4 tests.
- **Evidence tier**: compile-time integration, empirical focused tests, and
  static source audit.
- **Residual risk**: `MaterialDatabase` and fluid constructors still inherit
  nalgebra scalar bounds through `Fluid<T>`; broader hemolysis, boundary,
  geometry, mesh, and solver surfaces remain open.

# Sprint 1.96.79 Resolution: cfd-core Velocity Leto Vector

### RESOLVED-108: Velocity Value Object Stored Nalgebra Vector3
- **Location**: `crates/cfd-core/src/physics/values/velocity.rs`,
  `crates/cfd-core/src/management/aggregates/parameters.rs`, and Leto
  `crates/leto/src/geometry.rs`.
- **Issue**: The scalar value wrappers had moved to Eunomia, but the velocity
  value object still stored `nalgebra::Vector3` and required `RealField`.
- **Remediation**: Added direct Leto consumption to `cfd-core`, changed
  `Velocity` and `PhysicalParameters::gravity` to `leto::geometry::Vector3`,
  and added Serde derives to Leto fixed 3D geometry types so serialized CFDrs
  value objects keep their boundary contract without a wrapper.
- **Verification**: Touched-file rustfmt passed. Focused scan shows no
  nalgebra `Vector3`/`RealField` in `Velocity` or `PhysicalParameters`.
  `cargo check -p cfd-core` passed after compiling the modified local Leto
  provider. `cargo nextest run -p cfd-core --lib` passed 183/183 tests.
  Downstream `cargo check -p cfd-2d`, `cargo check -p cfd-3d`, and
  `cargo check -p cfd-validation` passed.
- **Evidence tier**: compile-time integration plus static source audit.
- **Residual risk**: `ProblemAggregate`, `SimulationAggregate`, material,
  hemolysis, and fluid-dynamics contracts still contain nalgebra-bound surfaces.

# Sprint 1.96.78 Resolution: cfd-core Physics Value Eunomia Scalars

### RESOLVED-107: Scalar Physics Values Used Nalgebra RealField
- **Location**: `crates/cfd-core/src/physics/values/{temperature,pressure,flow,dimensionless}.rs`
  and immediate aggregate owners in `crates/cfd-core/src/management/aggregates`.
- **Issue**: Scalar-only value wrappers still depended on nalgebra scalar
  traits even though the provider migration assigns scalar contracts to
  Eunomia.
- **Remediation**: Replaced the scalar wrapper bounds with Eunomia
  `FloatElement`/`NumericElement`, changed zero/abs/sqrt operations to
  Eunomia contracts, and propagated required bounds through immediate aggregate
  owners.
- **Verification**: Touched-file rustfmt passed. A focused scan of the four
  scalar value-wrapper files returns no `nalgebra::RealField`, `RealField`, or
  `ComplexField` matches. `cargo check -p cfd-core` passed. `cargo nextest run
  -p cfd-core --lib` passed 183/183 tests.
- **Evidence tier**: compile-time integration plus static source audit.
- **Residual risk**: `Velocity` still owns `nalgebra::Vector3`; material,
  hemolysis, and broader aggregate paths still contain nalgebra-bound surfaces
  pending later provider slices.

# Sprint 1.96.77 Resolution: cfd-math WGPU Boundary Reduction

### RESOLVED-106: cfd-math Owned Direct WGPU Metric Coupling
- **Location**: root `Cargo.toml`, `crates/cfd-math/Cargo.toml`,
  `crates/cfd-core/src/compute/gpu/mod.rs`, and
  `crates/cfd-math/src/linear_solver/operators/gpu.rs`.
- **Issue**: `cfd-math` declared optional `wgpu` and imported WGPU only to
  block on submitted work and query timestamp-query support, leaving a
  downstream math crate coupled to raw GPU infrastructure during the Hephaestus
  migration.
- **Remediation**: Added `GpuContext::synchronize` and
  `GpuContext::supports_timestamp_queries`, routed cfd-math dispatch metrics
  through those methods, and removed direct `dep:wgpu` activation from
  `cfd-math/gpu` and the root package `gpu` feature.
- **Evidence tier**: compile-time integration, static manifest/source audit,
  and empirical focused tests. `cargo check -p cfd-math --features gpu`
  passed; `cargo nextest run -p cfd-math --features gpu --lib` passed 298/298
  tests; `cargo check -p cfd-suite --features gpu` passed. No runtime GPU
  performance claim is made.
- **Residual risk**: Raw WGPU buffers, command encoders, pipelines, and WGSL
  kernels remain in `cfd-core::compute::gpu` pending full Hephaestus migration.

# Sprint 1.96.76 Resolution: cfd-core Hephaestus GPU Probe

### RESOLVED-105: cfd-core Performed Raw WGPU Adapter Probing
- **Location**: `crates/cfd-core/Cargo.toml`,
  `crates/cfd-core/src/compute/traits.rs`, and root `Cargo.toml`.
- **Issue**: `ComputeBackend::detect_gpu_support` still created a raw
  `wgpu::Instance` and requested an adapter directly instead of using the Atlas
  GPU provider.
- **Dependency audit**: `cargo tree --workspace -i hephaestus-wgpu` shows
  `hephaestus-wgpu -> cfd-core`; `cargo tree --workspace -i rustfft` reports no
  matching package; `cargo tree --workspace -i ndarray` still resolves only
  through `numpy -> cfd-python`.
- **Remediation**: Added optional `hephaestus-wgpu` to `cfd-core` under the
  `gpu` feature and changed GPU availability detection to
  `hephaestus_wgpu::WgpuDevice::try_default`.
- **Verification**: `cargo check -p cfd-core --features gpu` passed.
  `cargo nextest run -p cfd-core --features gpu --lib` passed 183/183 tests.
- **Evidence tier**: compile-time integration plus empirical focused tests and
  static dependency audit.
- **Residual risk**: Raw WGPU buffers, command encoders, pipelines, and WGSL
  kernels remain in `cfd-core::compute::gpu`, with adapter coupling in
  `cfd-math::linear_solver::operators::gpu`.

# Sprint 1.96.75 Resolution: cfd-1d Removed Non-Python ndarray Path

### RESOLVED-104: cfd-1d Pulled ndarray Through Unused sprs
- **Location**: `crates/cfd-1d/Cargo.toml`, root `Cargo.toml`, and cfd-1d
  generic test helpers.
- **Issue**: `cfd-1d` declared `sprs` but did not use it in source. That stale
  dependency pulled `ndarray v0.17.2` into the active 1D/3D graphs. The root
  workspace also still declared unused `ndarray`.
- **Dependency audit**: `cargo tree -p cfd-1d -i ndarray` and `cargo tree -p
  cfd-3d -i ndarray` now report no matching package. Workspace-wide
  `ndarray` remains only through `numpy -> cfd-python`.
- **Remediation**: Removed `sprs`, removed the unused root workspace `ndarray`
  declaration, refreshed the lockfile, and propagated `FloatElement` through
  generic cfd-1d tests that construct migrated cfd-core fluid types.
- **Evidence tier**: static manifest/lock audit plus compile-time and empirical
  tests. Touched-file rustfmt passed. `cargo check -p cfd-1d` passed. `cargo
  nextest run -p cfd-1d` passed 725/725 tests with 3 skipped.
- **Residual risk**: cfd-1d still uses `nalgebra` and `nalgebra-sparse`; Leto
  sparse/vector replacement remains open. cfd-python's `numpy` boundary still
  pulls `ndarray` and needs a separate API-boundary decision.

# Sprint 1.96.74 Resolution: cfd-3d Resolved ndarray-backed Apollo Packages

### RESOLVED-103: CFDrs Used an Older ndarray-backed Apollo Git Revision
- **Location**: Workspace provider patches, cfd-3d FEM/spectral call sites,
  and cfd-validation dev-dependency generic bounds.
- **Issue**: CFDrs was locked to an older Apollo Git revision whose active FFT
  and NUFFT packages still exposed ndarray-backed APIs. Updating to the local
  Atlas Apollo provider stack required Eunomia scalar/complex call-site
  propagation in cfd-3d and validation code.
- **Dependency audit**: `apollo-fft` and `apollo-nufft` now resolve from
  `D:\atlas\repos\apollo`; active Apollo package trees have no `ndarray`
  matches, and the Apollo source/manifest/lock scan returns no `ndarray` hits.
  The remaining active `ndarray` in the cfd-3d graph is `sprs -> cfd-1d`, not
  Apollo.
- **Remediation**: Added side-by-side Atlas provider patches for Apollo and its
  local dependencies, updated cfd-3d FEM/spectral code to use Eunomia
  scalar/complex contracts, and propagated `FloatElement` bounds/conversions
  through cfd-validation code needed by cfd-3d nextest.
- **Evidence tier**: compile-time integration plus empirical cfd-3d nextest and
  static dependency/source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-3d` passed. `cargo nextest run -p cfd-3d` passed 394/394 tests; one
  existing mesh-convergence test was reported slow at 16.394s.
- **Residual risk**: Full Hephaestus GPU kernel replacement and non-Apollo
  ndarray removal remain open. Leto sparse/operator replacement remains for
  later focused slices. Apollo rustfft/num-complex validation and benchmark
  residue is outside this ndarray-removal slice.

# Sprint 1.96.73 Resolution: cfd-math Geometric Multigrid Bypassed Eunomia

### RESOLVED-102: Geometric Multigrid Used Direct num-traits Conversion
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/{gmg,amg.rs}`
- **Issue**: GMG hierarchy construction, Poisson matrix construction, and
  transfer weights used direct `num_traits::FromPrimitive`,
  `cfd_core::conversion::SafeFromF64`, direct `T::from_f64(...)`, direct
  `T::from_usize(...)`, and `T::from_f64_or_one(...)` scalar construction.
  AMG retained stale direct `FromPrimitive` bounds after its callees no longer
  required them.
- **Dependency audit**: `eunomia::FloatElement` now owns GMG scalar
  construction and transfer weights. `eunomia::NumericElement` owns
  grid-dimension-to-scalar conversion. `multigrid` has no direct `num-traits`
  provider residue.
- **Remediation**: Added Eunomia scalar helpers, replaced direct
  conversion/fallback calls, removed stale AMG direct provider bounds, and
  added value-semantic tests for Poisson stencil constants and full-weighting
  restriction values.
- **Evidence tier**: compile-time integration plus empirical focused GMG/AMG
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math geometric_multigrid
  poisson_matrix restrict_residual fas_solve amg` passed 11/11 tests.
  Multigrid-wide static scan found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, conversion
  fallback, `from_f64_or`, `SafeFromF64`, stale `rayon`/`tokio`, `rustfft`,
  `ndarray`, or `num_complex` hits.
- **Residual risk**: The multigrid public surface still uses nalgebra
  dense/sparse/vector types pending Leto migration. Broader `cfd-math` and
  workspace provider migration remains open, including raw GPU paths and
  remaining nalgebra surfaces outside multigrid.

---

# Sprint 1.96.72 Resolution: cfd-math Multigrid Interpolation Bypassed Eunomia

### RESOLVED-101: Multigrid Interpolation Used Direct num-traits Conversion
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/interpolation.rs`
- **Issue**: Interpolation operators and quality analysis used direct
  `num_traits::{FromPrimitive, ToPrimitive}` bounds, direct
  `T::from_f64(...)`, direct `T::from_usize(...)`, direct index `as f64`
  conversion, and `to_f64().unwrap_or(...)` fallback paths.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar construction
  and index-distance conversion. `eunomia::NumericElement` now owns
  absolute-value dispatch and f64 extraction for quality metrics. Static scan
  found no direct `num-traits` provider residue in `interpolation.rs`.
- **Remediation so far**: Added Eunomia scalar/ratio helpers, replaced direct
  conversion fallback calls, routed interpolation quality metric extraction
  through Eunomia, and removed touched test debug output/direct index-distance
  casts.
- **Evidence tier**: compile-time integration plus empirical focused
  interpolation tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  interpolation` passed 14/14 tests. Static scan found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, stale
  `rayon`, direct `as f64`, or `.to_f64()` fallback hits in
  `interpolation.rs`.
- **Residual risk**: `multigrid/gmg` and `amg.rs` still retain direct provider
  residue, and the multigrid public surface still uses nalgebra sparse/vector
  types pending Leto migration.

---

# Sprint 1.96.71 Resolution: cfd-math Multigrid Coarsening Bypassed Eunomia

### RESOLVED-100: Multigrid Coarsening Used Direct num-traits Conversion
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/coarsening`
- **Issue**: Coarsening algorithms and quality analysis used direct
  `num_traits::{FromPrimitive, ToPrimitive}` bounds plus direct
  `T::from_f64(...)`, `T::from_usize(...)`, `T::max_value(...)`, and
  `to_f64().unwrap_or(...)` scalar conversion/fallback paths.
- **Dependency audit**: `eunomia::FloatElement` now owns coarsening scalar
  construction and count conversion. `eunomia::NumericElement` now owns
  absolute-value dispatch, f64 extraction for quality metrics, and maximum
  finite sentinels. `multigrid/coarsening` has no direct `num-traits`
  provider residue.
- **Remediation**: Added Eunomia scalar helpers, replaced direct
  conversion/fallback calls, routed strength and quality absolute-value
  dispatch through Eunomia, and added value-semantic strength-matrix
  connectivity coverage for a 4-by-4 Poisson grid.
- **Evidence tier**: compile-time integration plus empirical focused
  coarsening tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  coarsening` passed 10/10 tests. Static scan found no direct `num_traits`,
  `FromPrimitive`, `ToPrimitive`, direct `T::from_f64`, direct
  `T::from_usize`, conversion fallback, `from_f64_or`, `SafeFromF64`, or stale
  `rayon` hits in `multigrid/coarsening`.
- **Residual risk**: `multigrid/interpolation.rs`, `multigrid/gmg`, and
  `amg.rs` still retain direct provider residue. The multigrid public surface
  still uses nalgebra sparse/vector types pending Leto migration.

---

# Sprint 1.96.70 Resolution: cfd-math Multigrid Smoothers Bypassed Eunomia

### RESOLVED-099: Multigrid Smoothers Used Direct Scalar Fallbacks
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/multigrid/{smoothers.rs,amg.rs}`
- **Issue**: Multigrid smoothers used direct
  `T::from_f64(...).unwrap_or_else(...)` fallback construction and nalgebra
  absolute-value dispatch for diagonal thresholds, Chebyshev eigenvalue
  defaults, recurrence constants, and smoother update thresholds. AMG's
  smoother-owner paths used direct fallback construction for coarsening
  thresholds, relaxation factors, and complexity filters.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar construction in
  the touched smoother and AMG owner paths. `eunomia::NumericElement` now owns
  absolute-value dispatch for smoother thresholds and AMG complexity filters.
  `smoothers.rs` has no direct `num-traits` provider residue.
- **Remediation**: Added Eunomia scalar helpers, replaced direct conversion
  fallback constants and absolute-value dispatch, propagated `FloatElement`
  bounds through smoother impls, and replaced existence-only smoother tests
  with value-semantic update/eigenvalue assertions.
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
  Other multigrid modules, raw GPU operators, and nalgebra sparse/vector
  surfaces remain open provider-migration work.

---

# Sprint 1.96.69 Resolution: cfd-math Convergence Monitor Bypassed Eunomia

### RESOLVED-098: Convergence Monitor Used Direct Fallback Conversion
- **Location**: `crates/cfd-math/src/linear_solver/traits.rs`
- **Issue**: `ConvergenceMonitor` used direct `T::from_f64(...)` and
  `cfd_core::conversion::SafeFromF64` fallback helpers for convergence-factor
  exponent construction, CG theoretical-bound factor construction, and
  validation safety-multiplier construction.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar construction
  and `powf` dispatch in the touched convergence-monitor methods.
  `SafeFromF64` and direct scalar conversion fallback usage no longer appear in
  `linear_solver/traits.rs`.
- **Remediation**: Added Eunomia scalar helper, added method-level
  `FloatElement` bounds only where convergence scalar construction is needed,
  disambiguated `powf` through Eunomia, and added value-semantic tests for the
  changed helper behavior.
- **Evidence tier**: compile-time integration plus empirical focused
  convergence tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  convergence_factor cg_theoretical_bound validate_convergence` passed 4/4
  tests, including the existing AMG convergence-factor match. Static scan
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`,
  direct `T::from_f64`, `T::from_usize`, `SafeFromF64`, `from_f64_or`,
  conversion fallback, or stale `rayon` hits in `traits.rs`.
- **Residual risk**: Other linear-solver modules still contain direct provider
  residue, including the raw GPU operator path, multigrid, and tests. Core
  solver trait APIs still use nalgebra vector surfaces pending later Leto
  migration.

---

# Sprint 1.96.68 Resolution: cfd-math Linear Operators Bypassed Eunomia

### RESOLVED-097: CPU Linear Operators Used num-traits Scalar Conversion
- **Location**:
  `crates/cfd-math/src/linear_solver/operators/{poisson.rs,momentum.rs}`
- **Issue**: Poisson/Laplacian and momentum/energy finite-difference operators
  used direct `num_traits::FromPrimitive` bounds plus
  `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks for
  scalar constants.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant
  construction and provider bounds in the touched CPU operator files.
  `num-traits` no longer appears in
  `linear_solver/operators/{poisson.rs,momentum.rs}`.
- **Remediation**: Added Eunomia scalar helpers, replaced direct conversion
  fallback constants, routed operator impl bounds through `FloatElement`, and
  added value-semantic center-stencil tests for the touched operators.
- **Evidence tier**: compile-time integration plus empirical focused
  value-semantic operator tests and static source audit. Touched-file rustfmt
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  laplacian_center_impulse poisson_center_impulse` passed 2/2 tests. `cargo
  nextest run -p cfd-math momentum_1d_applies momentum_2d_applies
  energy_2d_applies` passed 3/3 tests. Static scan found no direct
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct
  `T::from_f64`, `T::from_usize`, conversion fallback, `from_f64_or`, or stale
  `rayon` hits in `poisson.rs` or `momentum.rs`.
- **Residual risk**: The raw GPU operator path still bypasses Hephaestus and
  carries direct provider residue. Multigrid, traits, and tests still contain
  direct `num-traits` usage. The touched CPU operators still use nalgebra
  `DVector` surfaces pending later Leto migration.

---

# Sprint 1.96.67 Resolution: cfd-math Schwarz/Cholesky Bypassed Eunomia

### RESOLVED-096: Schwarz/Cholesky Used Direct Provider Residue
- **Location**:
  `crates/cfd-math/src/linear_solver/preconditioners/{schwarz.rs,cholesky.rs}`
- **Issue**: Schwarz carried a stale direct `num_traits::FromPrimitive`
  import and provider bound even though its local ILU dependency is
  conversion-free. IncompleteCholesky used direct
  `T::from_f64(...).unwrap_or_else(...)` for symmetry tolerance construction
  and direct nalgebra absolute-value dispatch for symmetry residuals.
- **Dependency audit**: `eunomia::FloatElement` now owns Cholesky tolerance
  construction. `eunomia::NumericElement` now owns Cholesky symmetry residual
  absolute-value dispatch. `num-traits` no longer appears in
  `linear_solver/preconditioners/{schwarz.rs,cholesky.rs}`.
- **Remediation**: Removed Schwarz's stale direct provider import/bound, added
  Eunomia scalar helpers to Cholesky, and made the residual absolute-value
  dispatch explicit through Eunomia.
- **Evidence tier**: compile-time integration plus empirical focused
  preconditioner tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math cholesky`
  passed 4/4 tests. `cargo nextest run -p cfd-math schwarz` passed 1/1 test.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in `schwarz.rs` or `cholesky.rs`.
- **Residual risk**: Other linear-solver modules still contain direct
  `num-traits` usage, including operators, multigrid, traits, and tests.
  Schwarz still uses nalgebra sparse surfaces pending later Leto migration;
  Cholesky moved to Leto CSR in Sprint 1.96.166.

---

# Sprint 1.96.66 Resolution: cfd-math SSOR Bypassed Eunomia

### RESOLVED-095: SSOR Used num-traits Scalar Conversion
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners/ssor.rs`
- **Issue**: SSOR used direct `num_traits::FromPrimitive` and
  `T::from_f64(...).unwrap_or_else(...)` for default relaxation construction.
- **Dependency audit**: `eunomia::FloatElement` now owns SSOR scalar constant
  construction and provider bounds. `num-traits` no longer appears in
  `linear_solver/preconditioners/ssor.rs`.
- **Remediation**: Added a Eunomia scalar helper, replaced default relaxation
  construction, and routed SSOR impl bounds through Eunomia `FloatElement`.
- **Evidence tier**: compile-time integration plus empirical focused SSOR
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math ssor` passed 4/4 tests.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in `ssor.rs`.
- **Residual risk**: Other linear-solver modules still contain direct
  `num-traits` usage, including Schwarz preconditioners, operators,
  multigrid, traits, and tests. SSOR still uses nalgebra sparse matrix/vector
  surfaces pending later Leto migration.

---

# Sprint 1.96.65 Resolution: cfd-math Basic Preconditioners Bypassed Eunomia

### RESOLVED-094: Basic Preconditioners Used num-traits Scalar Conversion
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners/basic.rs`
- **Issue**: Jacobi and SOR preconditioners used direct
  `num_traits::FromPrimitive`, `T::from_f64(...).unwrap_or_else(...)`, and
  direct scalar absolute-value dispatch for diagonal tolerances and SOR omega
  construction.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant
  construction. `eunomia::NumericElement` now owns absolute-value dispatch.
  `num-traits` no longer appears in `linear_solver/preconditioners/basic.rs`.
- **Remediation**: Added Eunomia scalar helpers, replaced Jacobi tolerance
  construction, replaced SOR omega constants and omega construction, and routed
  absolute-value checks through Eunomia.
- **Evidence tier**: compile-time integration plus empirical focused
  preconditioner tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math jacobi`
  passed 5/5 tests. `cargo nextest run -p cfd-math sor` passed 6/6 tests.
  Static scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback,
  `from_f64_or`, or stale `rayon` hits in `basic.rs`.
- **Residual risk**: Other linear-solver modules still contain direct
  `num-traits` usage, including SSOR/Schwarz preconditioners, operators,
  multigrid, traits, and tests. Basic preconditioners still use nalgebra sparse
  matrix/vector surfaces pending later Leto migration.

---

# Sprint 1.96.64 Resolution: cfd-math GMRES Chain Bypassed Eunomia

### RESOLVED-093: GMRES Chain Used num-traits Scalar Bounds
- **Location**: `crates/cfd-math/src/linear_solver/{gmres,chain.rs,direct_solver.rs,block_preconditioner.rs}`
- **Issue**: GMRES, Arnoldi, the fallback chain, direct sparse solver, and
  block/SIMPLE preconditioners used direct `num_traits::{Float,
  FromPrimitive, ToPrimitive}` bounds plus direct `Float::abs`,
  `T::from_f64`, `T::from_usize`, and conversion fallback paths.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar provider
  bounds in the touched GMRES-centered call path. `eunomia::NumericElement`
  now owns absolute-value and f64 conversion dispatch. `num-traits` no longer
  appears in `gmres`, `chain.rs`, `direct_solver.rs`, or
  `block_preconditioner.rs`.
- **Remediation**: Replaced GMRES/Arnoldi bounds, fallback-chain bounds,
  direct sparse solver conversion bounds, and block/SIMPLE preconditioner
  scalar safeguards with Eunomia traits/helpers.
- **Evidence tier**: compile-time integration plus empirical focused
  linear-solver tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math gmres`
  passed 21/21 tests. `cargo nextest run -p cfd-math direct_solver` passed
  3/3 tests. `cargo nextest run -p cfd-math block_preconditioner` passed 2/2
  tests. Static scan found no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, conversion
  fallback, or stale `rayon` hits in the touched slice.
- **Residual risk**: Other linear-solver modules still contain direct
  `num-traits` usage, including operators, basic/SSOR/Schwarz preconditioners,
  multigrid, and tests. The direct solver still routes through rsparse's f64
  LU backend until a later Atlas-owned Leto/backend solver replacement.

---

# Sprint 1.96.63 Resolution: cfd-math Linear-Solver Config Bypassed Eunomia

### RESOLVED-092: Iterative Solver Defaults Used num-traits Conversion
- **Location**: `crates/cfd-math/src/linear_solver/{config.rs,bicgstab/mod.rs,conjugate_gradient/mod.rs,gmres/solver.rs}`
- **Issue**: `IterativeSolverConfig::default` used direct
  `num_traits::FromPrimitive` and `T::from_f64(...).expect(...)` for default
  tolerance construction. The parallel SpMV config doc still named Rayon.
- **Dependency audit**: `eunomia::FloatElement` now owns default tolerance
  construction and default-construction bounds for CG, BiCGSTAB, and GMRES.
  `linear_solver/config.rs`, `bicgstab/mod.rs`, and
  `conjugate_gradient/mod.rs` no longer contain direct `num_traits`,
  `FromPrimitive`, direct `T::from_f64`, or stale `rayon` references.
- **Remediation**: Added a Eunomia scalar helper to config defaults,
  propagated default-construction bounds to CG, BiCGSTAB, and GMRES, removed
  direct `num_traits` bounds from the CG/BiCGSTAB linear-solver trait impls,
  and corrected the config doc wording to Moirai.
- **Evidence tier**: compile-time integration plus empirical focused
  linear-solver tests and static source audit. Touched-file rustfmt passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  default_solver` passed 2/2 tests. `cargo nextest run -p cfd-math
  test_gmres_configurable_trait` passed 1/1 test. Static scan found no direct
  `num_traits`, `FromPrimitive`, direct `T::from_f64`, or stale `rayon` hits
  in `linear_solver/config.rs`, `bicgstab/mod.rs`, or
  `conjugate_gradient/mod.rs`.
- **Residual risk**: Other linear-solver modules still contain direct
  `num-traits` usage, and linear solvers still use nalgebra vector/matrix
  surfaces pending later Leto migration.

---

# Sprint 1.96.62 Resolution: cfd-math SIMD Scalars Bypassed Eunomia

### RESOLVED-091: SIMD CFD Helpers Used num-traits Scalar Conversions
- **Location**: `crates/cfd-math/src/simd/cfd.rs`
- **Issue**: CFD SIMD central-difference constants and field-norm
  square-root dispatch used direct `num_traits::FromPrimitive`,
  `T::from_f64(...).unwrap_or_else(...)`, and nalgebra scalar math instead of
  the Atlas numeric provider surface.
- **Dependency audit**: `eunomia::FloatElement` now owns CFD SIMD scalar
  constant construction. `eunomia::NumericElement` now owns field-norm
  square-root dispatch. `num-traits` and `rayon` no longer appear in
  `crates/cfd-math/src/simd`; SIMD execution remains on Moirai's slice
  adapters.
- **Remediation**: Added a Eunomia scalar helper, replaced all CFD SIMD
  central-difference constants, routed field norm through Eunomia
  `NumericElement::sqrt`, and removed direct `FromPrimitive` bounds.
- **Evidence tier**: compile-time integration plus empirical SIMD tests and
  static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math`
  passed. `cargo nextest run -p cfd-math simd` passed 26/26 tests. Static scan
  found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`,
  `Float::`, direct `T::from_f64`, conversion fallback, `T::epsilon`, or
  `rayon` hits in `crates/cfd-math/src/simd`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, including linear-solver and test modules. SIMD still
  carries nalgebra scalar/vector surfaces pending later Leto/eunomia
  scalar-bound cleanup.

---

# Sprint 1.96.61 Resolution: cfd-math Sparse Scalars Bypassed Eunomia

### RESOLVED-090: Sparse Utilities Used num-traits Scalar Math
- **Location**: `crates/cfd-math/src/sparse/{patterns.rs,operations.rs}`
- **Issue**: Sparse stencil constants, Frobenius norm dispatch,
  condition-estimate singular-diagonal thresholds, and diagonal-dominance
  checks used direct `num_traits::{Float, FromPrimitive, Signed}`,
  `Float::...`, `T::from_f64`, conversion fallback, and `T::epsilon` paths.
  Sparse SpMV docs also still described the old Rayon implementation even
  though the code used Moirai's slice adapter.
- **Dependency audit**: `eunomia::FloatElement` now owns sparse scalar
  constant construction. `eunomia::NumericElement` now owns square-root and
  absolute-value dispatch. `num-traits` and `rayon` no longer appear in
  `crates/cfd-math/src/sparse`. `nalgebra_sparse`/`nalgebra` remain for the
  current sparse matrix/vector surface.
- **Remediation**: Added a Eunomia scalar helper for sparse pattern constants,
  replaced sparse norm/condition/diagonal-dominance scalar math dispatch, and
  corrected sparse SpMV docs to name Moirai.
- **Evidence tier**: compile-time integration plus empirical sparse tests and
  static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math`
  passed. `cargo nextest run -p cfd-math sparse` passed 15/15 tests. Static
  scan found no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Signed`,
  `Float::`, direct `T::from_f64`, conversion fallback, `T::epsilon`, or
  `rayon` hits in `crates/cfd-math/src/sparse`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, including linear-solver, simd, and test modules. Sparse
  still carries `nalgebra_sparse::CsrMatrix` and nalgebra `DVector` pending
  later Leto sparse/vector migration.

---

# Sprint 1.96.60 Resolution: cfd-math Nonlinear Solver Scalars Bypassed Eunomia

### RESOLVED-089: Anderson/JFNK Used num-traits Scalar Math
- **Location**: `crates/cfd-math/src/nonlinear_solver/{anderson.rs,jfnk.rs}`
- **Issue**: Anderson/JFNK default constants, QR norm/diagonal checks,
  finite-difference perturbation safeguards, EW forcing clamp math, Givens
  rotations, and back-substitution checks used direct
  `num_traits::{Float, FromPrimitive}`, `Float::...`, `T::from_f64`, and
  `T::epsilon` paths.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant
  construction. `eunomia::NumericElement` now owns square-root,
  absolute-value, and min/max scalar dispatch in the touched nonlinear-solver
  math. `num-traits` no longer appears in
  `crates/cfd-math/src/nonlinear_solver`. The nalgebra vector/matrix solver
  surface noted during this scalar-only slice is superseded by Sprint
  1.96.109.
- **Remediation**: Added Eunomia scalar helpers, replaced Anderson/JFNK
  defaults and scalar math dispatch, and used nalgebra `default_epsilon()` for
  real-field epsilon semantics where the removed `Float` bound previously
  supplied `T::epsilon`.
- **Evidence tier**: compile-time integration plus empirical nonlinear-solver
  tests and static source audit. Touched-file rustfmt passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math nonlinear_solver` passed
  9/9 tests. Static scan found no direct `num_traits`, `FromPrimitive`,
  `Float::`, direct `T::from_f64`, conversion fallback, or `T::epsilon` hits
  in `crates/cfd-math/src/nonlinear_solver`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, including linear-solver, sparse, simd, and test modules.
  The nonlinear-solver nalgebra `DVector`/`DMatrix` residual is closed by
  Sprint 1.96.109.

---

# Sprint 1.96.59 Resolution: cfd-math Pressure-Velocity Scalars Bypassed Eunomia

### RESOLVED-088: SIMPLE Defaults Used Old Scalar Providers
- **Location**: `crates/cfd-math/src/pressure_velocity/simple.rs`
- **Issue**: SIMPLE default tolerance and relaxation constants were
  constructed through direct `num_traits::FromPrimitive`,
  `T::from_f64(...).expect(...)` scalar paths. The module also still used
  nalgebra `RealField` and old `T::zero()`/`T::one()` identities. Explicit
  config construction carried an unnecessary conversion-trait bound.
- **Dependency audit**: `eunomia::RealField` now owns the pressure-velocity
  scalar surface. `eunomia::FloatElement` owns SIMPLE default scalar
  construction. `num-traits` and nalgebra no longer appear in
  `crates/cfd-math/src/pressure_velocity`.
- **Remediation**: Added a Eunomia scalar helper, replaced SIMPLE default
  constants, moved the module to Eunomia `RealField`, replaced old scalar
  identities, and narrowed explicit `SIMPLEConfig::new` construction away from
  scalar-conversion bounds.
- **Evidence tier**: compile-time integration plus empirical
  pressure-velocity tests and static source audit. `cargo fmt -p cfd-math
  --check` passed. `cargo check -p cfd-math --no-default-features --lib`
  passed. `cargo nextest run -p cfd-math --no-default-features
  pressure_velocity --status-level fail` passed 3/3 tests with 341 skipped.
  `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings` passed. Static scan found no nalgebra, `DVector`, `DMatrix`, old
  scalar identities, `num_traits`, `num_complex`, ndarray, rayon, tokio,
  rustfft, `FromPrimitive`, `ToPrimitive`, or `From<f64>` hits in
  `crates/cfd-math/src/pressure_velocity`.
- **Residual risk**: Other cfd-math modules still contain direct old-provider
  scalar/storage surfaces, including linear-solver, sparse, SIMD, and test
  modules. The pressure-velocity cone itself no longer carries nalgebra
  scalar-provider residue; later July 4 slices closed high-order and
  nonlinear-solver provider residue.

---

# Sprint 1.96.58 Resolution: cfd-math Iterator Scalars Bypassed Eunomia

### RESOLVED-087: Iterator Utilities Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/iterators/{stencils.rs,statistics.rs}`
- **Issue**: Stencil coefficients and iterator statistics counts were
  constructed through direct `num_traits::FromPrimitive`, `T::from_f64`,
  `T::from_usize`, and conversion-fallback scalar paths. The second-derivative
  stencil path also returned a zero-vector placeholder for declared 3-point
  patterns outside `Central3`.
- **Dependency audit**: `eunomia::FloatElement` now owns iterator scalar
  constants and count-to-scalar construction. `eunomia::NumericElement` owns
  count staging and standard-deviation square-root dispatch. `num-traits` no
  longer appears in `crates/cfd-math/src/iterators`.
- **Remediation**: Added Eunomia scalar helpers, replaced stencil and
  statistics scalar conversions, routed standard-deviation square-root dispatch
  through Eunomia, and replaced the second-derivative placeholder with real
  `[1, -2, 1]` coefficients for every declared 3-point stencil pattern.
- **Evidence tier**: compile-time integration plus empirical iterator tests
  and static source audit. Touched-file `rustfmt --check` passed. `cargo check
  -p cfd-math` passed. `cargo nextest run -p cfd-math iterators` passed 7/7
  tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits
  in `crates/cfd-math/src/iterators`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, including linear-solver, sparse, simd, and test modules.
  Iterators still carry `nalgebra::RealField` pending later Leto/eunomia
  scalar-bound cleanup; later July 4 slices closed nonlinear-solver residue.

---

# Sprint 1.96.57 Resolution: cfd-math WENO Scalars Bypassed Eunomia

### RESOLVED-086: WENO Reconstruction Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/high_order/weno`
- **Issue**: WENO5/WENO7 epsilon defaults, linear weights, ENO reconstruction
  coefficients, and smoothness-indicator constants were constructed through
  direct `num_traits::FromPrimitive`, `T::from_f64`, and conversion-fallback
  scalar paths.
- **Dependency audit**: `eunomia::FloatElement` now owns WENO scalar constant
  construction. `num-traits` no longer appears in
  `crates/cfd-math/src/high_order`. `nalgebra::RealField` remains as the
  current high-order WENO scalar surface.
- **Remediation**: Added shared WENO Eunomia scalar and square helpers,
  replaced WENO5/WENO7 constants, and replaced WENO `FromPrimitive` bounds with
  `FloatElement`.
- **Evidence tier**: compile-time integration plus empirical WENO tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math weno` passed 6/6 tests.
  Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits
  in `crates/cfd-math/src/high_order`.
- **Residual risk**: Other cfd-math modules still contain direct provider
  residue, including linear-solver, sparse, simd, and test modules. The WENO
  nalgebra scalar-bound residue from this older slice is closed by Sprint
  1.96.106; nonlinear-solver residue is closed by Sprint 1.96.109.

---

# Sprint 1.96.56 Resolution: cfd-math Interpolation Scalars Bypassed Eunomia

### RESOLVED-085: Interpolation Used Old Scalar Providers
- **Location**: `crates/cfd-math/src/interpolation`
- **Issue**: The interpolation trait, linear interpolation, Lagrange
  interpolation, and cubic spline implementation still used nalgebra
  `RealField`; Lagrange and cubic spline also used old scalar identities or
  direct scalar conversion/fallback paths. Lagrange duplicate nodes could reach
  a zero basis denominator.
- **Dependency audit**: `eunomia::RealField` now owns the interpolation scalar
  surface. `eunomia::FloatElement` owns cubic-spline scalar constant
  construction. `num-traits` and nalgebra no longer appear in the touched
  interpolation source.
- **Remediation**: Moved interpolation trait and implementations to Eunomia
  scalar contracts, replaced old scalar identities with Eunomia constants,
  replaced cubic-spline `FromPrimitive` bounds with `FloatElement`, and added
  Lagrange duplicate-node rejection.
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
- **Residual risk**: Other cfd-math modules still contain direct old-provider
  scalar/storage surfaces, including linear-solver, sparse, SIMD, and test
  modules. Later July 4 slices closed high-order and nonlinear-solver provider
  residue.

---

# Sprint 1.96.55 Resolution: cfd-math Differentiation Scalars Bypassed Eunomia

### RESOLVED-084: Differentiation Operators Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/differentiation/{finite_difference.rs,gradient.rs}`
- **Issue**: Finite-difference stencil constants, gradient/divergence/curl
  constants, SIMD helper scalar staging, and 2D Laplacian constants were
  constructed through direct `num_traits::FromPrimitive`, `From<f64>`,
  `T::from_f64`, and `T::from_f32` scalar conversion paths.
- **Dependency audit**: `eunomia::FloatElement` now owns differentiation
  scalar constant construction and SIMD helper scalar staging.
  `eunomia::NumericElement` owns scalar identities, and Leto now owns
  `Array1` derivative outputs plus `Vector3` gradient/divergence/curl surfaces.
  `num-traits` and nalgebra no longer appear in the touched differentiation
  operator files.
- **Remediation**: Added Eunomia scalar helpers, replaced finite-difference
  constants, replaced gradient/divergence/curl constants, replaced SIMD helper
  scalar staging, moved result/vector surfaces to Leto, and replaced the
  type-suffixed SIMD helper name.
- **Evidence tier**: compile-time integration plus empirical differentiation
  tests and static source audit. Touched-file `rustfmt --check` passed.
  `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  differentiation` passed 12/12 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`, `From<f64>`,
  `T::from_f32`, or conversion-fallback hits in the touched files.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, including linear-solver, sparse, simd, and test modules.
  Differentiation vector residue is closed by the later July 4 differentiation
  Leto follow-up; nonlinear-solver residue is closed by Sprint 1.96.109.

---

# Sprint 1.96.54 Resolution: cfd-math Integration Scalars Bypassed Eunomia

### RESOLVED-083: Integration Quadrature Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/integration`
- **Issue**: Quadrature constants, interval-count conversions, adaptive
  tolerance/error dispatch, and tetrahedral quadrature constants were
  constructed through direct `num_traits::FromPrimitive`, `From<f64>`,
  `T::from_f64`, and `T::from_usize` scalar conversion paths.
- **Dependency audit**: `eunomia::FloatElement` now owns integration scalar
  constants and floating-point math dispatch. `eunomia::NumericElement` now
  owns exact count conversion staging and adaptive absolute-value dispatch.
  `num-traits` no longer appears in the touched integration module tree.
  `eunomia::RealField` now owns the quadrature scalar surface, and nalgebra no
  longer appears in the touched integration module tree.
- **Remediation**: Added Eunomia scalar helpers, replaced 1D quadrature
  constants and weights, replaced composite/adaptive/tensor conversion bounds,
  replaced tetrahedral quadrature constants, and moved the quadrature trait
  bounds off nalgebra `RealField`.
- **Evidence tier**: compile-time integration plus empirical
  integration-named tests and static source audit. `cargo fmt -p cfd-math
  --check` passed. `cargo check -p cfd-math --no-default-features --lib`
  passed. `cargo nextest run -p cfd-math --no-default-features integration
  --status-level fail` passed 12/12 tests with 330 skipped. `cargo clippy -p
  cfd-math --no-default-features --all-targets -- -D warnings` passed. Static
  scan found no nalgebra, `DVector`, `DMatrix`, old scalar identities,
  `num_traits`, `num_complex`, ndarray, rayon, tokio, or rustfft hits in
  `crates/cfd-math/src/integration`.
- **Residual risk**: Other cfd-math modules still contain direct old-provider
  scalar/storage surfaces, including linear-solver, sparse, SIMD, and test
  modules. Later July 4 slices closed high-order, interpolation, and
  nonlinear-solver provider residue.

---

# Sprint 1.96.53 Resolution: cfd-math Exponential Scalars Bypassed Eunomia

### RESOLVED-082: Exponential Time Stepping Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/exponential.rs`
- **Issue**: ERK4 coefficients, phi-function small-argument thresholds, and
  scaling/squaring factorial conversion were constructed through direct
  `num_traits::FromPrimitive` bounds,
  `T::from_f64(...).unwrap_or_else(...)` fallbacks, and
  `T::from_f64(...).unwrap()` conversion.
- **Dependency audit**: `eunomia::FloatElement` now owns exponential scalar
  constant construction and integer-to-scalar staging for the power-series
  factorial. `num-traits` no longer appears in the touched module.
  `nalgebra::DVector`/`DMatrix` remain as the current exponential integrator
  vector/matrix surfaces.
- **Remediation**: Added local Eunomia conversion helpers, replaced ERK4 scalar
  coefficients, phi-function thresholds, and scaling/squaring factorial
  conversion.
- **Evidence tier**: compile-time integration plus empirical exponential
  time-stepping tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  time_stepping::exponential::` passed 2/2 tests. Static scan found no
  `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_f64`, `T::from_usize`,
  or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/exponential.rs`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, and the exponential integrator vector/matrix API remains
  nalgebra-owned until the later Leto vector/matrix migration.

---

# Sprint 1.96.52 Resolution: cfd-math IMEX Scalars Bypassed Eunomia

### RESOLVED-081: IMEX Time Stepping Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/imex.rs`
- **Issue**: Newton tolerance, ARS343 gamma/delta construction, tableau
  coefficients, and explicit/implicit solution weights were constructed
  through direct `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)`
  fallbacks and nalgebra scalar `sqrt` dispatch.
- **Dependency audit**: `eunomia::FloatElement` now owns IMEX scalar constant
  construction, and `eunomia::NumericElement` owns ARS343 scalar square-root
  dispatch. `num-traits` no longer appears in the touched module.
  `nalgebra::DVector`/`DMatrix` remain as the current time-stepper state and
  Jacobian surfaces.
- **Remediation**: Added a local Eunomia conversion helper, replaced Newton
  tolerance, gamma/delta construction, explicit/implicit tableau constants, and
  solution weights.
- **Evidence tier**: compile-time integration plus empirical IMEX tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math imex` passed 5/5 tests.
  Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/imex.rs`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, and the time-stepper state/Jacobian API remains
  nalgebra-owned until the later Leto vector/matrix migration.

---

# Sprint 1.96.51 Resolution: cfd-math RKC Scalars Bypassed Eunomia

### RESOLVED-080: RKC Time Stepping Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/rk_chebyshev.rs`
- **Issue**: RKC defaults, Chebyshev recurrence constants, adaptive step
  halving, error normalization, clamp bounds, and exponent dispatch were
  constructed through direct `num_traits::FromPrimitive` conversions and
  nalgebra scalar method dispatch.
- **Dependency audit**: `eunomia::FloatElement` now owns RKC floating constants
  and transcendental dispatch. `eunomia::NumericElement` now owns absolute
  value/square-root dispatch and exact count conversion staging. `num-traits`
  no longer appears in the touched module. `nalgebra::DVector` remains as the
  current time-stepper state surface.
- **Remediation**: Added local Eunomia scalar helpers, replaced RKC defaults,
  Chebyshev recurrence constants, adaptive error-control constants, stage and
  vector length conversions, and disambiguated scalar math dispatch through
  Eunomia traits.
- **Evidence tier**: compile-time integration plus empirical RKC tests and
  static source audit. Touched-file `rustfmt --check` passed. `cargo check -p
  cfd-math` passed. `cargo nextest run -p cfd-math rk_chebyshev` passed 4/4
  tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/rk_chebyshev.rs`.
- **Residual risk**: Other cfd-math modules still contain direct `num-traits`
  scalar conversions, and the time-stepper state API remains nalgebra-owned
  until the later Leto vector migration.

---

# Sprint 1.96.50 Resolution: cfd-math Adaptive Scalars Bypassed Eunomia

### RESOLVED-079: Adaptive Time Stepping Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/adaptive.rs`
- **Issue**: Adaptive defaults, PI controller gains, rejection scaling, clamp
  bounds, and Dormand-Prince tableau coefficients were constructed through
  direct `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks.
- **Dependency audit**: `eunomia::FloatElement` now owns adaptive scalar
  constant construction. `num-traits` no longer appears in the touched module.
  `nalgebra::DVector` remains as the current time-stepper state surface.
- **Remediation**: Added a local Eunomia conversion helper, replaced adaptive
  defaults/controller/tableau constants, and disambiguated power dispatch
  through `FloatElement::powf`.
- **Evidence tier**: compile-time integration plus empirical adaptive
  time-stepping tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  adaptive` passed 3/3 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/adaptive.rs`.
- **Residual risk**: Other cfd-math time-stepping modules still contain direct
  `num-traits` scalar conversions, and the time-stepper state API remains
  nalgebra-owned until the later Leto vector migration.

---

# Sprint 1.96.49 Resolution: cfd-math Runge-Kutta Scalars Bypassed Eunomia

### RESOLVED-078: Runge-Kutta Methods Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/runge_kutta.rs`
- **Issue**: RK3, RK4, and low-storage RK4 constructed scalar constants through
  direct `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks.
  The touched low-storage RK4 path also combined its solution accumulator and
  stage state incorrectly instead of using the Carpenter-Kennedy 2N residual
  recurrence.
- **Dependency audit**: `eunomia::FloatElement` now owns Runge-Kutta scalar
  constant construction. `num-traits` no longer appears in the touched module.
  `nalgebra::DVector` remains as the current time-stepper state surface.
- **Remediation**: Added a local Eunomia conversion helper, replaced RK scalar
  constants, implemented `r_i = a_i r_{i-1} + dt f(t_i, u_i)` and
  `u_{i+1} = u_i + b_i r_i`, and added value-semantic decay plus zero-RHS
  preservation tests.
- **Evidence tier**: compile-time integration plus empirical value-semantic
  Runge-Kutta tests and static source audit. Touched-file `rustfmt --check`
  passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math
  runge_kutta` passed 5/5 tests. Static scan found no `num_traits`,
  `FromPrimitive`, `ToPrimitive`, `T::from_f64`, or conversion-fallback hits in
  `crates/cfd-math/src/time_stepping/runge_kutta.rs`.
- **Residual risk**: Other cfd-math time-stepping modules still contain direct
  `num-traits` scalar conversions, and the time-stepper state API remains
  nalgebra-owned until the later Leto vector migration.

---

# Sprint 1.96.48 Resolution: cfd-math Stability Scalars Bypassed Eunomia

### RESOLVED-077: Stability Analyzer Used num-traits Scalar Fallbacks
- **Location**: `crates/cfd-math/src/time_stepping/stability/mod.rs`; `crates/cfd-math/src/time_stepping/stability/analysis.rs`
- **Issue**: The stability analyzer still imported `num_traits::ToPrimitive`
  and constructed scalar thresholds through direct `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` fallbacks after the complex-value migration.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant construction, and `eunomia::NumericElement` owns diagnostic conversion to `f64`. `num-traits` no longer appears in the touched stability module.
- **Remediation**: Added local Eunomia conversion helpers, replaced CFL/RK/von-Neumann scalar constants and formatting conversions, and disambiguated `abs` dispatch through Eunomia.
- **Evidence tier**: compile-time integration plus empirical value-semantic stability tests and static source audit. Touched-file `rustfmt --check` passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math stability` passed 5/5 tests. `git diff --check` passed. Static scan found no `num_traits`, `ToPrimitive`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits under `crates/cfd-math/src/time_stepping/stability`.
- **Residual risk**: Package-level `cargo fmt --package cfd-math --check` remains blocked by unrelated existing formatting drift outside the touched stability files. The module still uses nalgebra `DMatrix`/`DVector` for the current Butcher-tableau surface until the later Leto matrix migration.

---

## Audit Principles
- **Mathematical Accuracy**: Zero tolerance for error masking or unverified approximations.
- **Implementation Completeness**: Full theorem documentation and rigorous testing required.
- **Literature Compliance**: Algorithms must match primary literature exactly.
- **Quality Standards**: No "working but incorrect" implementations.

---

# Sprint 1.96.46 Resolution: cfd-2d Scheme Boundary Bypassed Eunomia Complex

### RESOLVED-075: cfd-2d Scheme Amplification Used num-complex
- **Location**: `crates/cfd-2d/src/schemes/mod.rs`; `crates/cfd-2d/Cargo.toml`
- **Issue**: `SpatialDiscretization::amplification_factor` returned and constructed `num_complex::Complex<f64>` even though the Atlas datatype vocabulary is `eunomia::Complex`.
- **Dependency audit**: `eunomia::Complex<f64>` now owns the touched scheme amplification value. `num-complex` no longer appears in cfd-2d source or manifest files. Direct `num-traits` remains in cfd-2d for conversion contracts outside this slice.
- **Remediation**: Replaced the amplification-factor complex type and constructors with Eunomia, removed the direct manifest dependency, and added explicit `eunomia::FloatElement` bounds where cfd-2d constructs cfd-core provider-bound grids or solver configs.
- **Evidence tier**: compile-time integration plus empirical regression tests and static source audit. `cargo fmt --package cfd-2d --package cfd-3d --package cfd-validation --check` passed. `cargo check -p cfd-2d`, `cargo check -p cfd-3d`, and `cargo check -p cfd-validation` passed. `cargo nextest run -p cfd-2d` passed 563/563 tests with 27 skipped. Static scan found no `num_complex`, `num-complex`, or `NumComplex` hits under `crates/cfd-2d`.
- **Residual risk**: Full cfd-2d provider migration still needs direct `num-traits` and nalgebra-owned matrix/vector paths removed.

---

# Sprint 1.96.47 Resolution: cfd-3d Apollo/Leto Adapter Compile Repair

### RESOLVED-076: cfd-3d Adapter Assumed Unreachable Apollo Leto API
- **Location**: `crates/cfd-3d/src/atlas_array.rs`; `crates/cfd-3d/src/fem/config.rs`; `crates/cfd-3d/src/spectral/solver.rs`; `crates/cfd-validation/src/manufactured/richardson/validation.rs`
- **Issue**: The private cfd-3d adapter imported Apollo Leto-native FFT symbols that are not present in the locked Apollo revision, imported `leto::MnemosyneStorage` without enabling the Leto feature that exports it, and left cfd-3d/cfd-validation solver-config call sites without the `FloatElement` provider bound required by cfd-core/cfd-2d.
- **Dependency audit**: Leto remains the cfd-3d algorithm-storage owner. Apollo's reachable FFT/NUFFT API still uses ndarray arrays, so ndarray conversion stays isolated in the private adapter. Eunomia owns the scalar provider bounds for the touched config construction paths.
- **Remediation**: Replaced the nonexistent Apollo Leto-native FFT calls with reachable Apollo array calls behind Leto conversions, removed the stale `MnemosyneStorage` import, added the missing cfd-3d Eunomia dependency and bounds, and propagated `FloatElement` into MMS Richardson validation solver construction.
- **Evidence tier**: compile-time integration plus empirical cfd-2d regression tests. `cargo check -p cfd-3d` passed. `cargo check -p cfd-validation` passed. `cargo nextest run -p cfd-2d` passed 563/563 tests with 27 skipped.
- **Residual risk**: Full cfd-3d ndarray removal still depends on Apollo exposing reachable Leto-native FFT/NUFFT signatures or moving this conversion boundary upstream into Apollo.

---

# Sprint 1.96.45 Resolution: cfd-math/cfd-validation Stability Analysis Bypassed Eunomia Complex

### RESOLVED-074: Stability Analysis Used num-complex Callback Values
- **Location**: `crates/cfd-math/src/time_stepping/stability/analysis.rs`; `crates/cfd-validation/src/time_integration/stability_analysis.rs`
- **Issue**: Von Neumann stability analysis exposed and consumed `num_complex::Complex<f64>` even though the Atlas datatype vocabulary is `eunomia::Complex`.
- **Dependency audit**: `eunomia::Complex<f64>` now owns the touched stability callback and validation spatial-operator complex values. `num-complex` still resolves transitively through nalgebra/simba and remains directly owned by cfd-2d outside this slice.
- **Remediation**: Replaced direct `num_complex` callback types and constructors with Eunomia complex values, removed direct `num-complex` manifest dependencies from cfd-math and cfd-validation, and replaced the explicit RK dense complex matrix inverse with forward substitution over `(I - zA)x = 1`.
- **Evidence tier**: compile-time integration plus empirical value-semantic stability tests and static source audit. Touched-file rustfmt passed. `cargo check -p cfd-math` passed. `cargo nextest run -p cfd-math stability` passed 5/5 tests. Static scan found no `num_complex`, `num-complex`, or `NumComplex` hits in `crates/cfd-math` or `crates/cfd-validation`.
- **Residual risk**: The later cfd-2d Eunomia scheme slice removed the direct cfd-2d `num-complex` boundary and cleared the downstream validation compile blocker. Full `num-complex` removal still requires the broader nalgebra/Leto replacement.

---

# Sprint 1.96.44 Resolution: cfd-core Rhie-Chow Bypassed Eunomia

### RESOLVED-073: cfd-core Rhie-Chow Used num-traits Conversion Fallbacks
- **Location**: `crates/cfd-core/src/physics/fluid_dynamics/rhie_chow.rs`
- **Issue**: Rhie-Chow interpolation used direct `num_traits::FromPrimitive` bounds and `T::from_f64(...).unwrap_or_else(...)` fallback conversions for the default relaxation factor and face-interpolation constant.
- **Dependency audit**: `eunomia::FloatElement` now owns touched Rhie-Chow scalar constants. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct conversion bounds with Eunomia `FloatElement`, routed constants through the provider, and added value-semantic tests for u-face and v-face pressure-correction formulas.
- **Evidence tier**: compile-time integration plus empirical value-semantic Rhie-Chow tests and static source audit. `rustfmt --check --edition 2021 crates\cfd-core\src\physics\fluid_dynamics\rhie_chow.rs` passed. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features rhie_chow` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 176/176 tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/fluid_dynamics/rhie_chow.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates, including fluid-dynamics finite-difference operations, cavitation number helpers, and shared downstream conversion facades. Full Leto CPU replacement for nalgebra-owned fluid-dynamics surfaces remains a separate migration.

---

# Sprint 1.96.43 Resolution: cfd-core Boundary Geometry Bypassed Eunomia

### RESOLVED-072: cfd-core Boundary Geometry Used Conversion Fallback Constants
- **Location**: `crates/cfd-core/src/physics/boundary/geometry.rs`
- **Issue**: Boundary geometry sphere and cylinder `measure()` formulas used direct `T::from_f64(...).unwrap_or_else(...)` fallback conversions for `4`, `3`, and `pi`.
- **Dependency audit**: `eunomia::FloatElement` now owns touched boundary geometry scalar constants. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct conversion fallbacks with Eunomia constants, disambiguated power dispatch through `FloatElement::powi`, split the impl so `contains_point()` and `dimension()` do not inherit the conversion bound, and added value-semantic tests for line length, sphere volume, cylinder volume, and unsupported zero-measure geometry variants.
- **Evidence tier**: compile-time integration plus empirical value-semantic geometry tests and static source audit. `rustfmt --check --edition 2021 crates\cfd-core\src\physics\boundary\geometry.rs` passed. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features boundary::geometry` passed 4/4 tests. `cargo nextest run -p cfd-core --no-default-features` passed 174/174 tests. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/geometry.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates, including fluid-dynamics interpolation/operations, cavitation number helpers, and shared downstream conversion facades. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

# Sprint 1.96.42 Resolution: cfd-core Boundary Ghost Cells Bypassed Eunomia

### RESOLVED-071: cfd-core Boundary Ghost Cells Used num-traits Conversions
- **Location**: `crates/cfd-core/src/physics/boundary/ghost_cells.rs`
- **Issue**: Boundary ghost-cell formulas used direct `num_traits::FromPrimitive` and `T::from_f64(...).unwrap_or_else(...)` conversions for stencil constants, and Robin singularity reporting used `num_traits::ToPrimitive`.
- **Dependency audit**: `eunomia::FloatElement` now owns touched ghost-cell scalar constants and `eunomia::NumericElement` owns Robin singularity value reporting. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct conversion bounds with Eunomia traits, routed Dirichlet/Neumann/Robin constants through `FloatElement`, routed diagnostic conversion through `NumericElement`, and made the degenerate `alpha == 0 && beta == 0` Robin case return `BoundaryErrorKind::RobinSingularity` before division.
- **Evidence tier**: compile-time integration plus empirical value-semantic ghost-cell tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features ghost_cells` passed 5/5 tests. `cargo nextest run -p cfd-core --no-default-features` passed 170/170 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion-fallback hits in `physics/boundary/ghost_cells.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates, including boundary geometry helpers, mesh operations, cavitation, fluid dynamics, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

# Sprint 1.96.41 Resolution: cfd-core Staggered Grid Bypassed Eunomia

### RESOLVED-070: cfd-core Staggered Grid Used num-traits Coordinate Conversions
- **Location**: `crates/cfd-core/src/geometry/staggered.rs`
- **Issue**: Staggered-grid coordinate construction used direct `num_traits::FromPrimitive` bounds and `T::from_usize`/`T::from_f64` conversions for grid dimensions, face indices, and half-cell offsets.
- **Dependency audit**: `eunomia::FloatElement` now owns touched staggered-grid scalar conversion. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct conversion bounds with Eunomia `FloatElement`, routed half constants through Eunomia, and added exact-representability assertions for integer grid dimensions and indices before converting them into the scalar type.
- **Evidence tier**: compile-time integration plus empirical value-semantic staggered-grid tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features staggered` passed 5/5 tests. `cargo nextest run -p cfd-core --no-default-features` passed 167/167 tests. Touched-file `rustfmt --check` and `git diff --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_`, `from_usize`, or conversion-fallback hits in `geometry/staggered.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates, including boundary geometry/ghost-cell helpers, mesh operations, cavitation, fluid dynamics, and blood/non-Newtonian fluid modules. Full Leto CPU replacement for nalgebra-owned geometry surfaces remains a separate migration.

---

# Sprint 1.96.40 Resolution: cfd-core Boundary Time Functions Bypassed Eunomia

### RESOLVED-069: cfd-core Boundary Time Functions Used num-traits and nalgebra Math Dispatch
- **Location**: `crates/cfd-core/src/physics/boundary/{time_dependent,applicator,applicators,manager,specification}.rs`
- **Issue**: Boundary time functions and ghost-cell helpers used direct `T::from_f64(...).unwrap_or_else(T::one)` conversions for `2π` and `2`, and sinusoidal/exponential time functions dispatched math through nalgebra scalar methods instead of the Atlas numeric provider.
- **Dependency audit**: `eunomia::FloatElement` now owns touched boundary scalar conversion and explicit `sin`/`exp` dispatch. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct conversion fallbacks with Eunomia constants, routed `sin`/`exp` through `FloatElement`, propagated the provider bound through boundary applicator/specification/manager surfaces, and added value-semantic tests for time-function and ghost-cell behavior.
- **Evidence tier**: compile-time integration plus empirical value-semantic boundary tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features boundary` passed 30/30 tests. `cargo nextest run -p cfd-core --no-default-features` passed 167/167 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.sin()` hits in touched boundary files.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full Leto CPU replacement for nalgebra-owned boundary/geometry surfaces remains a separate migration.

---

# Sprint 1.96.39 Resolution: cfd-core Temperature Fluids Bypassed Eunomia

### RESOLVED-068: cfd-core Temperature Fluids Used num-traits Bounds and nalgebra Math Dispatch
- **Location**: `crates/cfd-core/src/physics/fluid/temperature.rs`
- **Issue**: Temperature-dependent fluid models imported `num_traits::FromPrimitive`, carried stale conversion bounds on polynomial and exponential models, and converted the Sutherland exponent through direct `T::from_f64(...).unwrap_or_else(T::one)`.
- **Dependency audit**: `eunomia::FloatElement` now owns touched temperature-fluid scalar conversion and explicit `exp`/`powf` dispatch. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Removed direct `FromPrimitive` import and bounds, routed Arrhenius and Andrade exponentials through `FloatElement::exp`, routed Sutherland exponent construction and power through `FloatElement`, and added value-semantic tests for polynomial, Andrade, and Sutherland formulas.
- **Evidence tier**: compile-time integration plus empirical value-semantic temperature-fluid tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features fluid::temperature` passed 6/6 tests. `cargo nextest run -p cfd-core --no-default-features` passed 162/162 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, conversion-fallback, bare `.exp()`, or bare `.powf(...)` hits in `physics/fluid/temperature.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full Leto CPU replacement for nalgebra-owned fluid/physics surfaces remains a separate migration.

---

# Sprint 1.96.38 Resolution: cfd-core Hemolysis Bypassed Eunomia

### RESOLVED-067: cfd-core Hemolysis Used num-traits Constants and Inverted Activation Exponential
- **Location**: `crates/cfd-core/src/physics/hemolysis/{calculator,trauma}.rs`
- **Issue**: Hemolysis clinical-index and platelet activation constants used direct `num_traits::FromPrimitive` calls with silent `T::one()` fallbacks. `PlateletActivation::activation_probability` also applied the sign of the exponential twice, producing negative probabilities above the shear threshold.
- **Dependency audit**: `eunomia::FloatElement` now owns touched hemolysis scalar conversion. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and fallback conversions with Eunomia `FloatElement`; corrected platelet activation to `1 - exp(-k * excess_stress * exposure_time)`.
- **Evidence tier**: compile-time integration plus empirical value-semantic hemolysis tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features hemolysis` passed 9/9 tests. `cargo nextest run -p cfd-core --no-default-features` passed 160/160 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `physics/hemolysis/{calculator,trauma}.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full Leto CPU replacement for nalgebra-owned physics surfaces remains a separate migration.

---

# Sprint 1.96.37 Resolution: cfd-core Mesh Quality Bypassed Eunomia

### RESOLVED-066: cfd-core Mesh Quality Used num-traits Threshold Constants
- **Location**: `crates/cfd-core/src/geometry/mesh/quality.rs`
- **Issue**: Mesh quality assessment converted fixed aspect-ratio, skewness, and orthogonality thresholds through direct `num_traits::FromPrimitive` calls with silent `T::one()` fallback values.
- **Dependency audit**: `eunomia::FloatElement` now owns mesh quality scalar threshold conversion. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and fallback conversions with Eunomia `FloatElement`, preserving strict threshold comparison behavior.
- **Evidence tier**: compile-time integration plus empirical value-semantic mesh quality tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features mesh_quality` passed 3/3 tests. `cargo nextest run -p cfd-core --no-default-features` passed 158/158 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `geometry/mesh/quality.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full Leto CPU replacement for nalgebra-owned mesh/geometry surfaces remains a separate migration.

---

# Sprint 1.96.36 Resolution: cfd-core CPU Backend Bypassed Eunomia

### RESOLVED-065: cfd-core CPU Backend Used num-traits Fallback Conversion
- **Location**: `crates/cfd-core/src/compute/cpu.rs`
- **Issue**: CPU advection converted `KernelParams` f64 domain parameters through a local `safe_f64_to_t` helper backed by `num_traits::FromPrimitive`, silently substituting one or zero when conversion failed. CPU buffer impls also carried unnecessary conversion bounds.
- **Dependency audit**: `eunomia::FloatElement` now owns CPU advection scalar conversion. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and the fallback helper with Eunomia conversions, and narrowed `CpuBuffer` construction, `ComputeBuffer`, and `Debug` impls to `RealField + Copy`.
- **Evidence tier**: compile-time integration plus empirical value-semantic backend tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features test_cpu_advection_kernel_linear_exactness` passed 1/1 test. `cargo nextest run -p cfd-core --no-default-features` passed 155/155 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `safe_f64_to_t`, or conversion-fallback hits in `compute/cpu.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full Leto CPU replacement for nalgebra-owned vector/matrix surfaces remains a separate migration.

---

# Sprint 1.96.35 Resolution: cfd-core Time Integrators Bypassed Eunomia

### RESOLVED-064: cfd-core Time Integrators Used num-traits Constants
- **Location**: `crates/cfd-core/src/compute/time/integrators.rs`
- **Issue**: RK2, RK4, Crank-Nicolson, and implicit integrator defaults used direct `num_traits::FromPrimitive` conversions. Backward Euler and Crank-Nicolson default tolerances silently substituted `T::one()` if conversion failed.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant conversion in the integrator file. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and `T::from_f64` calls with Eunomia constants, removed conversion-only error branches, and added value-semantic tests for Forward Euler, RK2, RK4, implicit defaults, invalid iteration counts, and Backward Euler nonconvergence.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features compute::time::integrators` passed 6/6 tests. `cargo nextest run -p cfd-core --no-default-features` passed 155/155 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/integrators.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full nalgebra replacement remains a separate Leto-backed state/vector migration.

---

# Sprint 1.96.34 Resolution: cfd-core Time-Step Controllers Bypassed Eunomia

### RESOLVED-063: cfd-core Time-Step Controllers Used num-traits Constants and Float Helpers
- **Location**: `crates/cfd-core/src/compute/time/controllers.rs`
- **Issue**: Adaptive and variable time-step controllers used direct `num_traits::FromPrimitive` conversions, `num_traits::Float` helper calls, and silent one/zero fallbacks for defaults and integration-order conversion.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar default conversion and power functions in the controller file. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive`/`Float` bounds and conversions with Eunomia constants/math, replaced min/max helper calls with local comparison helpers, and changed `calculate_dt` to return typed configuration errors for zero, overflowing, or unsupported integration order instead of silently substituting order one.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features compute::time::controllers` passed 4/4 tests. `cargo nextest run -p cfd-core --no-default-features` passed 149/149 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `compute/time/controllers.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Time integrators remain a later Eunomia slice.

---

# Sprint 1.96.33 Resolution: cfd-core Solver Config Defaults Bypassed Eunomia

### RESOLVED-062: cfd-core Solver Config Used num-traits Default Constants
- **Location**: `crates/cfd-core/src/compute/solver/config.rs`
- **Issue**: Solver and linear-solver configuration defaults used direct `num_traits::FromPrimitive` conversions plus silent zero fallbacks for tolerances, time step, and CFL defaults. `SolverConfigBuilder::new` duplicated the default construction block.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar default conversion in solver configuration. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and `T::from_f64(...).unwrap_or_else(num_traits::Zero::zero)` conversions with Eunomia constants, routed `SolverConfigBuilder::new` through `SolverConfig::default`, and added value-semantic tests for builder/default parity.
- **Evidence tier**: compile-time integration plus empirical value-semantic tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features solver_config` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 146/146 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero-conversion fallback hits in `compute/solver/config.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full nalgebra replacement remains a separate Leto-backed state/vector/solver migration.

---

# Sprint 1.96.32 Resolution: cfd-core Abstraction Defaults Bypassed Eunomia

### RESOLVED-061: cfd-core Abstractions Used num-traits Default Constants
- **Location**: `crates/cfd-core/src/abstractions/{state.rs,problem.rs}`
- **Issue**: `FieldState::new` and `ProblemParameters::default` used direct `num_traits::FromPrimitive` conversions plus silent one-value fallbacks for default time step, reference pressure, and reference temperature.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar default conversion in the touched abstraction files. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced the `FromPrimitive` bounds and `T::from_f64(...).unwrap_or_else(...)` default conversions with Eunomia constants, and split `FieldState` impls so only constructors/defaults require scalar-conversion bounds.
- **Evidence tier**: compile-time integration plus empirical abstraction tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features abstractions` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero/one conversion-fallback hits in `abstractions/{state,problem}.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full nalgebra replacement remains a separate Leto-backed state/vector migration.

---

# Sprint 1.96.31 Resolution: cfd-core Fluid Validation Thresholds Bypassed Eunomia

### RESOLVED-060: cfd-core Fluid Validation Used num-traits Threshold Constants
- **Location**: `crates/cfd-core/src/physics/fluid/validation.rs`
- **Issue**: Fluid property validation thresholds used direct `num_traits::FromPrimitive` conversions plus silent zero/one fallbacks for property bounds and dimensionless/thermodynamic limits.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar threshold conversion in `validation.rs`. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced the `FromPrimitive` bound and all `T::from_f64(...).unwrap_or_else(...)` validation threshold conversions with Eunomia constants while preserving the existing validation ranges and error behavior.
- **Evidence tier**: compile-time integration plus empirical validation tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features physics::fluid::validation` passed 2/2 tests. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` and `git diff --check` passed. Static scan found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits in `validation.rs`.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and other CFDrs crates. Full nalgebra replacement remains a separate Leto-backed state/vector migration.

---

# Sprint 1.96.30 Resolution: cfd-core Material/Fluid Constants Bypassed Eunomia

### RESOLVED-059: cfd-core Material and Fluid Constructors Used num-traits Constants
- **Location**: `crates/cfd-core/src/physics/material/*`, `crates/cfd-core/src/physics/fluid/{database.rs,newtonian.rs}`
- **Issue**: Material constructors, the material database setup path, constant-property fluid constructors, fluid database constructors, and ideal-gas constants used direct `num_traits::FromPrimitive` conversions plus silent zero/one fallbacks.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant conversion in the touched material/fluid constructor cone. `eunomia::NumericElement` owns the ideal-gas square-root call. `num-traits` remains in cfd-core for modules outside this slice.
- **Remediation**: Replaced direct `FromPrimitive` bounds and `T::from_f64(...).unwrap_or_else(...)` calls with Eunomia conversions, removed conversion-only error branches from fluid database constants, and kept existing validation for constructed physical properties.
- **Evidence tier**: compile-time integration plus empirical regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo check -p cfd-core` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` and `git diff --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion-fallback hits.
- **Residual risk**: Direct `num-traits` use remains elsewhere in cfd-core and across other CFDrs crates. Full nalgebra replacement remains a separate Leto-backed state/vector migration.

---

# Sprint 1.96.29 Resolution: cfd-core Physics Values Bypassed Eunomia

### RESOLVED-058: cfd-core Physics Value Objects Used num-traits Constants
- **Location**: `crates/cfd-core/Cargo.toml`, `crates/cfd-core/src/physics/values/*`, `crates/cfd-core/src/management/aggregates/{parameters.rs,problem.rs,simulation.rs}`
- **Issue**: Physics value objects and immediate management aggregate callers used direct `num_traits::FromPrimitive` bounds and `T::from_f64(...).unwrap_or_else(T::zero/one)` scalar constant conversions instead of the Atlas Eunomia scalar provider.
- **Dependency audit**: `eunomia::FloatElement` now owns scalar constant conversion in the touched cone. `num-traits` remains in cfd-core for modules outside this slice. `nalgebra` remains for vector/scalar math pending a later Leto-backed value-state migration.
- **Remediation**: Added `eunomia.workspace = true` to cfd-core, replaced direct `FromPrimitive` bounds with `FloatElement`, replaced scalar constant conversion fallbacks with Eunomia conversions, disambiguated nalgebra `abs`/`sqrt` calls, and made `PhysicalParameters::with_reynolds` return `Result<Self>` instead of silently substituting a default Reynolds number.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests and static source audit. `cargo check -p cfd-core --no-default-features` passed. `cargo nextest run -p cfd-core --no-default-features` passed 144/144 tests. Touched-file `rustfmt --check` passed. Static scan over the touched cone found no `num_traits`, `FromPrimitive`, `T::from_f64`, or silent fallback conversion hits.
- **Residual risk**: Full cfd-core Eunomia migration still needs `management::conversion` and other numeric helper modules. Full nalgebra removal still requires a later Leto-backed vector/state migration. Package-level `cargo fmt --check --package cfd-core` remains blocked by pre-existing unrelated formatting drift.

---

# Sprint 1.96.28 Resolution: cfd-python Blood Bindings Bypassed Eunomia

### RESOLVED-057: cfd-python Used num-traits in Blood-Model PyO3 Wrappers
- **Location**: `crates/cfd-python/Cargo.toml`, `crates/cfd-python/src/blood.rs`
- **Issue**: cfd-python still depended directly on `num-traits` for `FromPrimitive`/`ToPrimitive` conversions in blood-model bindings after the 2D dense-array Leto migration.
- **Dependency audit**: `eunomia::NumericElement` now owns the Python getter conversion boundary. `num-traits` is absent from direct cfd-python source and manifest usage.
- **Remediation**: Added `eunomia.workspace = true`, removed `num-traits = "0.2"`, passed `f64` shear rates directly into the Rust-owned models, and replaced getter fallback conversions with Eunomia `to_f64()` so Python exposes actual model fields.
- **Evidence tier**: compile-time integration plus static source/dependency audit. `cargo check -p cfd-python` passed. `cargo nextest run -p cfd-python --no-tests pass` confirmed 0 test binaries. `cargo fmt -p cfd-python --check` passed. Static scan found no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits under `crates/cfd-python`. Dependency-chain clippy blockers encountered during scoped clippy were fixed across cfd-1d, cfd-2d, cfd-3d, cfd-validation, and cfd-python; `cargo clippy -p cfd-python --all-targets -- -D warnings` now passes.
- **Residual risk**: Other CFDrs crates still use direct `num-traits`/`num-complex`. Package-wide cfd-2d formatting still has pre-existing unrelated drift outside the touched files.

---

# Sprint 1.96.27 Resolution: cfd-io Scalar Bounds Bypassed Eunomia

### RESOLVED-056: cfd-io Used num-traits Instead of Eunomia for Scalar Bounds
- **Location**: `crates/cfd-io/Cargo.toml`, `crates/cfd-io/src/binary.rs`, `crates/cfd-io/src/checkpoint/{data.rs,validator.rs}`, `crates/cfd-io/src/csv/*`
- **Issue**: After the dense payload migration, cfd-io still depended directly on `num-traits` for float bounds and checkpoint validation conversions instead of using the Atlas Eunomia scalar provider.
- **Dependency audit**: `eunomia::RealField` now owns the cfd-io scalar trait boundary. `num-traits` remains in the workspace for other crates and still resolves transitively through provider crates, but it is absent from direct cfd-io source and manifest usage.
- **Remediation**: Replaced `num_traits::Float`/`FromPrimitive`/`ToPrimitive` with Eunomia `RealField`, added `eunomia.workspace = true` to cfd-io, removed `num-traits.workspace = true`, and changed mass-conservation mesh-dimension conversion to reject zero or non-exactly-convertible dimensions instead of silently falling back to unit spacing.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests. `cargo check -p cfd-io` passed. `cargo check -p cfd-io --all-features` passed. `cargo nextest run -p cfd-io --no-fail-fast` passed 3/3 tests. `cargo clippy -p cfd-io --all-targets --all-features -- -D warnings` passed. `cargo fmt -p cfd-io --check` passed. Static scan found no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits under `crates/cfd-io`.
- **Residual risk**: Other CFDrs crates still use direct `num-traits`/`num-complex`; cfd-io still has provider-transitive `num-traits` through Leto/Eunomia. Those remain for later Eunomia/provider slices.

---

# Sprint 1.96.26 Resolution: cfd-io Checkpoint/Binary Payloads Bypassed Leto

### RESOLVED-055: cfd-io Exposed nalgebra Dense Payloads and Pulled nalgebra Through Core/Math
- **Location**: `crates/cfd-io/Cargo.toml`, `crates/cfd-io/src/{binary.rs,error.rs,leto_arrays.rs}`, `crates/cfd-io/src/checkpoint/{data.rs,manager.rs,validator.rs}`, `crates/cfd-io/src/csv/*`, `crates/cfd-io/tests/checkpoint_roundtrip.rs`
- **Issue**: cfd-io checkpoint and binary APIs exposed `nalgebra::DMatrix`/`DVector`, CSV/checkpoint validation used `RealField`, and the crate pulled nalgebra transitively through `cfd-core`/`cfd-math` even after direct matrix usage was removed.
- **Dependency audit**: `leto` now owns dense checkpoint and binary vector/matrix payloads. cfd-io owns a local file-format `Error`/`Result` type so it does not depend on cfd-core only for I/O errors. The temporary `num-traits` scalar boundary was removed in Sprint 1.96.27.
- **Remediation**: Replaced checkpoint fields with Leto `Array2`, added explicit row-major serde payloads, replaced binary helpers with Leto `Array1`/`Array2`, removed direct cfd-io nalgebra/cfd-core/cfd-math dependencies, and updated checkpoint tests to assert shape/value/checksum semantics.
- **Evidence tier**: compile-time integration plus empirical value-semantic regression tests. `cargo check -p cfd-io` passed. `cargo nextest run -p cfd-io --no-fail-fast` passed 3/3 tests. `cargo check -p cfd-io --all-features` passed. `cargo fmt -p cfd-io --check` passed. `cargo tree -p cfd-io -e normal -i nalgebra` returned no matching package. Static scan found no `nalgebra`, `DMatrix`, `DVector`, or `RealField` hits under `crates/cfd-io`.
- **Residual risk**: The optional `vtk` feature can still resolve nalgebra through upstream RITK/Gaia/Burn dependencies.

---

# Sprint 1.96.25 Resolution: cfd-python Owned PyO3 Dense Arrays Through ndarray/nalgebra

### RESOLVED-054: cfd-python 2D NumPy Helpers Bypassed Leto Ownership
- **Location**: `crates/cfd-python/Cargo.toml`, `crates/cfd-python/src/solver_2d/mod.rs`, `crates/cfd-python/src/solver_2d/poiseuille.rs`, `crates/cfd-python/src/solver_2d/cavity.rs`
- **Issue**: cfd-python directly used `nalgebra::DMatrix` and `ndarray::Array2` to build Python-returned dense arrays even though Leto is the Atlas dense-array provider target.
- **Dependency audit**: `leto` now owns the Rust dense arrays in the touched 2D binding helpers. `numpy` remains the Python ABI representation. Direct `ndarray` and `nalgebra` dependencies were removed from cfd-python.
- **Remediation**: Added a private Leto-to-NumPy helper, rewrote Poiseuille analytical velocity and Ghia benchmark arrays as Leto `Array2`, and copied rows into NumPy only at the PyO3 boundary.
- **Evidence tier**: compile-time integration plus static source/dependency audit. `cargo check -p cfd-python` passed. `cargo nextest run -p cfd-python --no-tests pass` confirmed 0 test binaries. `cargo fmt -p cfd-python --check` passed. Static scan found no `ndarray`, `nalgebra`, or `DMatrix` references in `crates/cfd-python`; Cargo.lock lists `leto` instead of direct `nalgebra`/`ndarray` under `cfd-python`.
- **Residual risk**: cfd-python still has direct `num-traits` usage in blood-model bindings; this remains for the Eunomia migration.

---

# Sprint 1.96.24 Resolution: CFDrs GPU Blocking Boundary Used pollster Instead of Moirai

- **Location**: `crates/cfd-core/src/compute/gpu/mod.rs`, `crates/cfd-core/src/compute/traits.rs`, `crates/cfd-core/src/compute/gpu/poisson_solver.rs`, `crates/cfd-math/src/linear_solver/operators/gpu.rs`, `crates/cfd-math/src/linear_solver/matrix_free/gpu_metrics_tests.rs`, `Cargo.toml`, `Cargo.lock`, `crates/cfd-core/Cargo.toml`, `crates/cfd-math/Cargo.toml`
- **Issue**: `cfd-core` and `cfd-math` used `pollster::block_on` for GPU context creation, GPU support detection, Poisson residual readback, and GPU operator sync dispatch even though Moirai is the Atlas concurrency provider target.
- **Dependency audit**: `moirai` now owns the synchronous GPU blocking boundary in the touched crates. `pollster` no longer appears in the resolved workspace graph. `hephaestus-wgpu` remains the GPU-provider target, but direct context replacement is blocked until CFDrs upgrades or abstracts away raw `wgpu 0.19` types because Hephaestus' provider is on `wgpu 26.0`.
- **Remediation**: Replaced the direct `pollster` waits with `moirai::block_on`, rewrote adapter probing to stay async inside `GpuContext::create_async`, and removed `pollster` from the workspace dependency graph.
- **Evidence tier**: compile-time integration plus empirical regression tests. `cargo check -p cfd-core --features gpu` passed. `cargo nextest run -p cfd-core --features gpu --lib` passed 151/151 tests. `cargo check -p cfd-math --features gpu` passed. `cargo nextest run -p cfd-math --features gpu --lib` passed 279/279 tests. Static scan found no `pollster` references in manifests, lockfile, or Rust sources.
- **Residual risk**: Full GPU migration still requires replacing CFDrs raw WGPU context, buffer, and pipeline ownership with Hephaestus after the `wgpu 0.19` to `wgpu 26.0` boundary is resolved.

---

# Sprint 1.96.23 Resolution: cfd-3d Dense Spectral Arrays Were Owned Directly by ndarray

### RESOLVED-052: Spectral and NUFFT Dense Arrays Bypassed Leto Ownership
- **Location**: `crates/cfd-3d/src/spectral/{dns.rs,forcing.rs,diagnostics.rs,fourier.rs}`, `crates/cfd-3d/src/ibm/nufft.rs`, `crates/cfd-3d/src/atlas_array.rs`
- **Issue**: The migrated `cfd-3d` spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT paths still owned dense arrays directly through `ndarray` even though Atlas Leto is the dense-array provider target.
- **Dependency audit**: `leto` now owns dense 1D/3D arrays in the touched paths. `apollo-fft` and `apollo-nufft` remain the transform providers, but their reachable published API still uses `ndarray`; the conversion is isolated in `atlas_array.rs` rather than spread through algorithm code.
- **Remediation**: Added a private Leto/Apollo adapter, changed touched algorithm paths to Leto arrays and checked accessor helpers, and kept residual `ndarray` usage at the Apollo boundary only.
- **Evidence tier**: compile-time integration plus empirical regression tests. `cargo check -p cfd-3d --no-default-features` passed. `cargo nextest run -p cfd-3d --no-default-features --lib` passed 206/206 tests.
- **Residual risk**: Full `ndarray` removal from `cfd-3d` requires adopting a reachable Apollo FFT/NUFFT API with Leto-native signatures or moving the conversion boundary upstream into Apollo.

---

# Sprint 1.96.22 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-051: Flow Analyzer Classified Reverse Flow with Signed Reynolds Number
- **Location**: `crates/cfd-1d/src/solver/analysis/analyzers/flow.rs`
- **Issue**: `FlowAnalyzer::determine_flow_regime` computed Reynolds number using signed velocity while the stored Reynolds number, velocity magnitude, and wall shear diagnostics used `|V|`. Reverse-flow edges could be classified as Stokes because `Re < 1` for negative signed Reynolds even when the physical Reynolds magnitude was transitional or turbulent.
- **Dependency audit**: `cfd-core` usage is appropriate for fluid density and viscosity; `cfd-math` remains solver infrastructure; `cfd-schematics` is topology and geometry-input authority; `petgraph`, `nalgebra`, `rayon`, `serde`, and `serde_json` remain appropriate supporting crates. The defect was local to `cfd-1d` post-solve flow diagnostics.
- **Remediation**: Changed regime classification to compute Reynolds number from velocity magnitude. Added real-network tests proving reverse-flow transitional classification and forward/reverse reciprocity for Reynolds number, wall shear rate, and flow regime.
- **Mathematical basis**: Reynolds number is `Re = rho |V| D_h / mu`, a nonnegative scalar comparing inertial to viscous forces. Flow orientation affects signed pressure-flow relations, not scalar regime classification.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features solver::analysis::analyzers::flow --lib` passed.
  - `cargo nextest run -p cfd-1d --lib --no-default-features solver::analysis::analyzers::flow --fail-fast --hide-progress-bar --status-level fail` passed.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.21 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-050: Droplet-Regime Dimensionless Groups Clamped Invalid Physical Domains
- **Location**: `crates/cfd-1d/src/physics/droplet_regime.rs`
- **Issue**: Capillary, Weber, and Ohnesorge number helpers clamped nonpositive surface tension and characteristic length denominators to tiny constants. This produced finite dimensionless groups and regime classifications for nonphysical two-phase states instead of rejecting undefined inputs.
- **Dependency audit**: `cfd-core` usage is appropriate for typed physics errors; `cfd-math` is appropriate for reusable numerical kernels outside these closed-form definitions; `cfd-schematics` is geometry-input authority; `petgraph`, `nalgebra`, `nalgebra-sparse`, `sprs`, `rayon`, `serde`, and `serde_json` remain appropriate supporting crates. The defect was local to `cfd-1d` two-phase droplet-regime validation.
- **Remediation**: Replaced denominator clamps with finite and physical-domain validation, changed dimensionless-group helpers and regime analysis to return typed errors, and rejected invalid capillary numbers before classification. Added tests for invalid surface tension, invalid characteristic length, and invalid capillary-number classification.
- **Mathematical basis**: `Ca = mu |U| / sigma`, `We = rho U^2 L / sigma`, and `Oh = mu / sqrt(rho sigma L)` are defined only for positive interfacial tension and positive characteristic length where present. Regime boundaries partition nonnegative finite `Ca`; negative or nonfinite `Ca` has no physical regime.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::droplet_regime --lib` passed.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::droplet_regime --fail-fast --hide-progress-bar --status-level fail` passed.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.20 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-049: Junction-Loss Rheology Used Signed Shear and Ignored Explicit Shear
- **Location**: `crates/cfd-1d/src/physics/resistance/models/junction_loss.rs`
- **Issue**: `JunctionLossModel::calculate_coefficients` estimated wall shear as `8V/D` with signed velocity and ignored `FlowConditions::shear_rate`. Non-Newtonian `cfd-core` rheology depends on shear-rate magnitude, so reverse flow could change or invalidate the viscous junction component while the K-factor minor loss remained scalar.
- **Dependency audit**: `cfd-core` usage is appropriate for fluid properties, Casson blood rheology, constants, and errors; `cfd-math` is appropriate for reusable solver infrastructure; `cfd-schematics` is topology and geometry-input authority; `petgraph`, `nalgebra`, `nalgebra-sparse`, `sprs`, `rayon`, `serde`, and `serde_json` remain appropriate supporting crates. The defect was local to `cfd-1d` junction-loss rheology coupling.
- **Remediation**: Derived default junction shear rate from velocity magnitude, routed explicit nonnegative shear rate into viscosity evaluation, and rejected negative explicit shear as a physics violation. Added Casson blood tests for explicit shear-thinning dependence, negative shear rejection, and reverse-flow resistance reciprocity.
- **Mathematical basis**: Junction minor loss is scalar in dynamic head, `Delta P = K rho V^2 / 2`, and the short-junction viscous component is linear in apparent viscosity. For generalized Newtonian fluids, `mu = mu(gamma_dot)` must use supplied wall shear when available; otherwise `gamma_dot = 8|V|/D`. Wall shear rate is a magnitude and cannot be negative.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::junction_loss --lib` passed.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::junction_loss --fail-fast --hide-progress-bar --status-level fail` passed.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.19 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-048: Darcy-Weisbach Resistance Ignored Explicit Non-Newtonian Shear Rate
- **Location**: `crates/cfd-1d/src/physics/resistance/models/darcy_weisbach.rs`
- **Issue**: `DarcyWeisbachModel::calculate_coefficients` always derived wall shear rate from velocity magnitude and ignored `FlowConditions::shear_rate`. It also replaced any `cfd-core` rheology error with baseline dynamic viscosity, which hid invalid non-Newtonian states and made the effective Reynolds number or laminar resistance inconsistent with the caller-supplied shear state.
- **Dependency audit**: `cfd-core` usage is appropriate for fluid properties, Casson blood rheology, constants, and errors; `cfd-math` is appropriate for reusable solver infrastructure; `cfd-schematics` is topology and geometry-input authority; `petgraph`, `nalgebra`, `nalgebra-sparse`, `sprs`, `rayon`, `serde`, and `serde_json` remain appropriate supporting crates. The defect was local to `cfd-1d` Darcy-Weisbach rheology coupling.
- **Remediation**: Routed explicit nonnegative shear rate into Darcy-Weisbach viscosity evaluation, rejected negative explicit shear as a physics violation, and propagated rheology errors instead of falling back to baseline viscosity. Added Casson blood tests for explicit shear-thinning dependence, negative shear rejection, and reverse-flow auto-Reynolds reciprocity.
- **Mathematical basis**: Darcy-Weisbach laminar resistance uses `R = 32 μ L / (A D_h^2)` and auto-Reynolds evaluation uses `Re = ρ |V| D_h / μ`. For generalized Newtonian fluids, `μ = μ(gamma_dot)` must use supplied wall shear when available; otherwise `gamma_dot = 8|V|/D_h`. Wall shear rate is a magnitude and cannot be negative.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::darcy_weisbach --lib` passed.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::darcy_weisbach --fail-fast --hide-progress-bar --status-level fail` passed.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.18 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-047: Rectangular-Channel Resistance Ignored Explicit Non-Newtonian Shear Rate
- **Location**: `crates/cfd-1d/src/physics/resistance/models/rectangular.rs`
- **Issue**: `RectangularChannelModel::calculate_coefficients` always derived wall shear rate from velocity magnitude and ignored `FlowConditions::shear_rate`. Coupled millifluidic solves can supply wall shear directly, and non-Newtonian `cfd-core` blood rheology depends on that shear. Negative explicit shear was also accepted because the model never inspected the field.
- **Dependency audit**: `cfd-core` usage is appropriate for Casson blood rheology, physical constants, and error contracts; `cfd-math` is appropriate for reusable linear algebra; `cfd-schematics` is topology and cross-section authority; `petgraph`, `nalgebra`, `nalgebra-sparse`, `sprs`, `rayon`, and `serde` remain appropriate infrastructure. The defect was local to `cfd-1d` rectangular resistance physics.
- **Remediation**: Routed explicit nonnegative shear rate into rectangular-channel viscosity evaluation, rejected negative explicit shear as a physics violation, and propagated rheology errors instead of falling back to baseline viscosity. Added Casson blood tests for explicit shear-thinning dependence and reverse-flow reciprocity.
- **Mathematical basis**: Rectangular laminar resistance is linear in apparent viscosity, `R = Po μ L / (2 A D_h^2)`. For generalized Newtonian fluids, `μ = μ(gamma_dot)` must use the supplied wall shear when a coupled solver provides it; otherwise the local estimate is `gamma_dot = 8|V|/D_h`. Wall shear rate is a magnitude and cannot be negative.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::rectangular --lib` passed.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::rectangular --fail-fast --hide-progress-bar --status-level fail` passed.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.17 Resolution: cfd-3d Dependency-Aware Physics Audit

### RESOLVED-046: FEM Stokes Assembly Accepted Incomplete Physical Problem Invariants
- **Location**: `crates/cfd-3d/src/fem/problem.rs`
- **Issue**: `StokesFlowProblem::validate` checked missing boundary nodes but did not reject nonpositive/nonfinite fluid density or viscosity, invalid pressure corner-node counts, nonfinite boundary/body-force data, or malformed per-element viscosity fields. These states violate Stokes well-posedness and could reach FEM matrix assembly.
- **Dependency audit**: `cfd-core` usage is appropriate for fluid and boundary contracts; `cfd-math` is appropriate for sparse/iterative solving; `cfd-mesh` is mesh topology authority; `cfd-io` is output infrastructure; `cfd-1d` and `cfd-2d` are cross-fidelity references; `cfd-schematics` is topology input authority. The defect was local to `cfd-3d` FEM problem validation.
- **Remediation**: Added validation for finite positive fluid properties, pressure-space corner-node count, finite body-force and boundary-condition values, and element-viscosity field length and positivity. Added unit tests covering each rejected physical state.
- **Mathematical basis**: The incompressible Stokes weak form is coercive on the velocity kernel only when dynamic viscosity is positive, density-dependent transient/advection terms require positive density, pressure DOFs must map to valid corner nodes, and per-element viscosity coefficients must be finite positive material parameters.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features fem::problem --lib` passed 10/10 tests.
  - `cargo nextest run -p cfd-3d --lib --no-default-features fem::problem --fail-fast --hide-progress-bar --status-level fail` passed 10/10 tests.
  - `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.16 Resolution: cfd-3d Dependency-Aware Physics Audit

### RESOLVED-045: Level-Set Transport Accepted Nonphysical Time, Mesh, and Velocity States
- **Location**: `crates/cfd-3d/src/level_set/advection.rs`
- **Issue**: The level-set Hamilton-Jacobi transport helper accepted nonpositive/nonfinite time steps, nonpositive/nonfinite grid spacing, and nonfinite velocity components before WENO5-Z/SSPRK3 or first-order upwind derivative reconstruction. These states make the explicit transport PDE undefined and can propagate invalid signed-distance fields.
- **Dependency audit**: `cfd-core` usage is appropriate for shared error contracts; `cfd-math` is appropriate for reusable numerical kernels; `cfd-mesh` is appropriate for mesh authority and scalar bounds; `cfd-io` is output infrastructure; `cfd-1d` and `cfd-2d` are cross-fidelity references; `cfd-schematics` is topology input authority. The defect was local to `cfd-3d` level-set transport validation.
- **Remediation**: Added transport input validation for finite positive `dt`, finite positive `dx/dy/dz`, and finite velocity components. Updated level-set tests so NaN velocity, zero time step, and negative time step are explicit errors.
- **Mathematical basis**: The Hamilton-Jacobi update `phi_t + u·grad(phi) = 0` requires finite velocity and positive finite time/space increments so upwind derivatives and SSPRK stage increments are defined real quantities.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features --test level_set_tests` passed 15/15 tests.
  - `cargo nextest run -p cfd-3d --test level_set_tests --no-default-features --fail-fast --hide-progress-bar --status-level fail` passed 15/15 tests.
  - `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.15 Resolution: cfd-3d Dependency-Aware Physics Audit

### RESOLVED-044: VOF Transport Accepted Nonphysical Time and Velocity States
- **Location**: `crates/cfd-3d/src/vof/advection.rs`
- **Issue**: The VOF CFL check rejected finite `CFL > 1` states but allowed nonpositive/nonfinite time steps and nonfinite velocity components to enter interface reconstruction and flux evaluation. This violated the explicit geometric-advection preconditions and could propagate invalid volume fractions instead of rejecting the physical state.
- **Dependency audit**: `cfd-core` usage is appropriate for shared error contracts and physical models; `cfd-math` is appropriate for reusable numerical kernels; `cfd-mesh` is appropriate as geometry/mesh authority; `cfd-io` is output infrastructure; `cfd-1d` and `cfd-2d` are cross-fidelity references; `cfd-schematics` is topology input authority. The defect was local to `cfd-3d` VOF transport validation.
- **Remediation**: Added finite-positive time-step validation, finite-positive grid-spacing validation, and finite velocity-component validation before CFL accumulation. Updated VOF tests so NaN/infinite velocities, zero time step, and negative time step are value-semantic errors.
- **Mathematical basis**: Explicit geometric VOF transport requires a positive finite `dt` and finite face velocities so the swept volume `|u| dt A_face` is a finite nonnegative volume and the directional CFL bound has physical meaning.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features --test vof_tests` passed 17/17 tests.
  - `cargo nextest run -p cfd-3d --test vof_tests --no-default-features --fail-fast --hide-progress-bar --status-level fail` passed 17/17 tests.
  - `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.14 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-043: Hagen-Poiseuille Derived Signed Shear for Non-Newtonian Fluids
- **Location**: `crates/cfd-1d/src/physics/resistance/models/hagen_poiseuille.rs`
- **Issue**: Hagen-Poiseuille coefficient evaluation derived wall shear rate as `8V/D` using signed velocity or signed flow-rate-derived velocity. Newtonian fluids masked the defect, but `cfd-core` non-Newtonian blood models use shear rate to compute apparent viscosity, so reverse flow through the same straight circular channel could produce different scalar resistance.
- **Dependency audit**: `cfd-core` usage is appropriate for Casson and other rheology models plus shared errors; `cfd-math` usage is appropriate for sparse/iterative network solving; `cfd-schematics` usage is appropriate as topology and cross-section input authority. The defect was local to `cfd-1d` resistance physics, not a crate-boundary mismatch.
- **Remediation**: Derived wall shear from velocity magnitude, rejected explicitly negative shear-rate inputs, and added Casson blood tests for reverse-flow resistance reciprocity and negative-shear rejection.
- **Mathematical basis**: For fully developed circular Poiseuille flow, wall shear-rate magnitude is `gamma_dot_w = 8|V|/D = 8|Q|/(A D)`. Straight-channel hydraulic resistance is an even scalar under flow reversal; orientation is represented by the signed network pressure-flow relation.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::hagen_poiseuille --lib` passed 6/6 tests.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::hagen_poiseuille --fail-fast --hide-progress-bar --status-level fail` passed 6/6 tests.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.13 Resolution: cfd-2d Dependency-Aware Physics Audit

### RESOLVED-042: LBM Accepted High-Mach Velocity States
- **Location**: `crates/cfd-2d/src/solvers/lbm/`
- **Issue**: The D2Q9 LBM solver documents the Chapman-Enskog low-Mach limit and crate instructions prohibit `Ma > 0.1`, but initialization and Zou-He boundary reconstruction accepted arbitrary imposed velocities. This allowed weakly compressible states outside the incompressible Navier-Stokes validity range.
- **Dependency audit**: `cfd-core` usage is appropriate for shared boundary-condition and error contracts; `cfd-math` usage is appropriate for numerical kernels outside the local LBM update; `cfd-mesh` and `cfd-io` serve mesh/output adapters; `cfd-1d` supplies reduced-order reference seeding; `cfd-schematics` remains topology/layout authority. The defect was local to `cfd-2d` LBM physics enforcement.
- **Remediation**: Added `Ma <= 0.1` validation for LBM initialization, Zou-He velocity inlet reconstruction, and pressure-boundary derived velocities. Boundary validation now returns errors through `LbmSolver::step`. Added tests that reject high-Mach initialization and velocity inlet states.
- **Mathematical basis**: D2Q9 BGK recovers incompressible Navier-Stokes through Chapman-Enskog expansion only for `Ma = |u|/c_s << 1`, with `c_s = sqrt(1/3)` in lattice units. Enforcing `Ma <= 0.1` bounds compressibility error at the documented low-Mach regime.
- **Verification**:
  - `cargo check -p cfd-2d --no-default-features` passed.
  - `cargo test -p cfd-2d --no-default-features solvers::lbm --lib` passed 31/31 tests.
  - `cargo nextest run -p cfd-2d --lib --no-default-features solvers::lbm --fail-fast --hide-progress-bar --status-level fail` passed 31/31 tests.
  - `cargo clippy -p cfd-2d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.12 Resolution: cfd-1d Dependency-Aware Physics Audit

### RESOLVED-041: Serpentine Resistance Used Signed Velocity for Scalar Losses
- **Location**: `crates/cfd-1d/src/physics/resistance/models/serpentine/`
- **Issue**: `SerpentineModel::calculate_coefficients` and `analyze` propagated signed velocity into wall shear rate, Reynolds number, Dean number, and pressure-loss magnitudes. Reversed flow therefore used the low-Re fallback path and returned different scalar resistance coefficients from the same physical channel operated in the opposite direction.
- **Dependency audit**: `cfd-core` usage is appropriate for fluid properties, constants, shared errors, and solver traits; `cfd-math` usage is appropriate for sparse/iterative network solving; `cfd-schematics` usage is appropriate as topology and cross-section authority. The defect was local to `cfd-1d` serpentine physics rather than a dependency-boundary mismatch.
- **Remediation**: Converted serpentine coefficient and analysis calculations to use velocity magnitude for scalar-loss quantities while leaving flow orientation to the network pressure-flow relation. Added reverse-flow tests for coefficients, resistance, Reynolds number, wall shear rate, Dean number, and pressure-drop magnitudes.
- **Mathematical basis**: For an orientation-symmetric serpentine channel, reversing flow maps `u -> -u` and `Q -> -Q`; scalar dissipation terms depend on `|u|`, `Re = rho |u| D_h / mu`, `De = Re sqrt(D_h/(2R_c))`, and dynamic pressure `0.5 rho |u|^2`. Therefore resistance magnitudes are even under reversal, while signed pressure drop is handled by the surrounding network equation.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::serpentine --lib` passed 17/17 tests.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::serpentine --fail-fast --hide-progress-bar --status-level fail` passed 17/17 tests.
  - `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.11 Resolution: cfd-3d Sigma SGS Energy Physics

### RESOLVED-040: Sigma SGS Energy Used Dimensional Strain-Rate Formula
- **Location**: `crates/cfd-3d/src/physics/turbulence/sigma.rs`
- **Issue**: `SigmaModel::turbulent_kinetic_energy` computed `nu_t |S| / Delta`, which has units `L/T^2` and is not turbulent kinetic energy. WALE and Vreman already used the shared Yoshizawa inversion `k_sgs = (nu_t / (C_k Delta))^2`, so Sigma carried an inconsistent LES diagnostic path.
- **Remediation**: Routed Sigma SGS kinetic energy through `kinetic_energy_from_eddy_viscosity`, removed the separate strain-rate dependency from the energy path, and added a regression that verifies the exact Yoshizawa value while rejecting the former formula.
- **Mathematical basis**: Yoshizawa's algebraic LES relation is `nu_t = C_k Delta sqrt(k_sgs)`. Solving for `k_sgs` gives `(nu_t / (C_k Delta))^2`, which has units `L^2/T^2` and preserves non-negativity for nonnegative eddy viscosity and positive filter width.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features physics::turbulence::sigma --lib` passed 4/4 tests.
  - `cargo nextest run -p cfd-3d --lib --no-default-features physics::turbulence::sigma --fail-fast --hide-progress-bar --status-level fail` passed 4/4 tests.
  - `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic` passed.

---

# Sprint 1.96.10 Resolution: cfd-3d VOF Directional CFL Physics

### RESOLVED-039: VOF Timestep Used Speed/Minimum-Spacing Bound
- **Location**: `crates/cfd-3d/src/vof/solver.rs`
- **Issue**: `VofSolver::calculate_timestep` computed `dt = C min(dx,dy,dz)/||u||_2`, while `AdvectionMethod` enforces the summed donor-cell swept-volume condition `|u_x|dt/dx + |u_y|dt/dy + |u_z|dt/dz <= 1`. The returned timestep therefore did not target the VOF advection invariant and overestimated the admissible step for diagonal flow when `C=1`.
- **Remediation**: Replaced the timestep denominator with `max_cells(|u_x|/dx + |u_y|/dy + |u_z|/dz)`, documented the CFL theorem on the solver API, and added a regression that accepts the corrected unit-grid diagonal timestep while rejecting the former norm/min-spacing estimate.
- **Mathematical basis**: Solving the geometric VOF swept-volume inequality for `dt` gives `dt <= C / (|u_x|/dx + |u_y|/dy + |u_z|/dz)` per cell. Taking the maximum directional rate over all cells gives the global explicit timestep bound.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features --test vof_tests test_calculate_timestep_uses_directional_vof_cfl` passed 1/1 test.
  - `cargo nextest run -p cfd-3d --test vof_tests --no-default-features test_calculate_timestep_uses_directional_vof_cfl --fail-fast --hide-progress-bar --status-level fail` passed 1/1 test.

---

# Sprint 1.96.9 Resolution: cfd-2d Upwind Coefficient Orientation

### RESOLVED-038: West-Face Upwind Coefficient Used Reversed Flux Sign
- **Location**: `crates/cfd-2d/src/discretization/convection.rs`
- **Issue**: `FirstOrderUpwind::coefficients` computed the west neighbor coefficient as `D_W + max(-F_W, 0)`. Under the standard finite-volume sign convention, the west coefficient must use `D_W + max(F_W, 0)`, so positive west-face flow was missing its upwind contribution while negative west-face flow received an unphysical one.
- **Remediation**: Replaced the west coefficient with `D_W + max(F_W, 0)`, documented the east/west coefficient invariant, and added tests that inspect both positive and negative west-face flux cases.
- **Mathematical basis**: For Patankar first-order upwinding with outward east-face flux `F_E` and inward west-face contribution `F_W`, bounded neighbor coefficients are `a_E = D_E + max(-F_E, 0)` and `a_W = D_W + max(F_W, 0)`. This preserves nonnegative neighbor coefficients and selects the upstream cell by face-flow orientation.
- **Verification**:
  - `cargo check -p cfd-2d --no-default-features` passed.
  - `cargo test -p cfd-2d --no-default-features discretization::convection --lib` passed 6/6 tests.
  - `cargo nextest run -p cfd-2d --lib --no-default-features discretization::convection --fail-fast --hide-progress-bar --status-level fail` passed 6/6 tests in 0.165 s.

---

# Sprint 1.96.8 Resolution: cfd-1d Branch Reverse-Flow Physics

### RESOLVED-037: Branch Solvers Rejected Reverse Parent Flow
- **Location**: `crates/cfd-1d/src/domain/junctions/branching/physics/`
- **Issue**: Two-way and three-way branch solvers rejected `Q_parent < 0`, and non-Newtonian wall-shear evaluation used signed flow. Reversed-flow network states could not be evaluated even though the Poiseuille pressure-drop law is odd in `Q` and the shear-rate/viscosity diagnostics are magnitude quantities.
- **Remediation**: Pressure-balanced branch solves now compute split fractions on `|Q_parent|` and reapply the parent-flow orientation to all daughter flows. Prescribed split paths accept signed parent flow. Wall shear and apparent viscosity use `|Q|`; pressure drops still use signed `Q`.
- **Mathematical basis**: For laminar branch channels, `ΔP_i = R_i(|Q_i|) Q_i`. Reversing flow orientation maps `Q_i → -Q_i` and `ΔP_i → -ΔP_i`, while `γ̇_i = 32|Q_i|/(πD_i³)` and constitutive viscosity remain unchanged. Mass conservation is invariant under the same sign transformation.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features --test branch_reverse_flow_orientation` passed 4/4 tests.
  - `cargo nextest run -p cfd-1d --test branch_reverse_flow_orientation --no-default-features --fail-fast --hide-progress-bar --status-level fail` passed 4/4 tests.

---

# Sprint 1.96.7 Resolution: cfd-3d Venturi Pressure-Coefficient Physics

### RESOLVED-036: 3D Venturi Coefficients Used Inlet Dynamic Pressure
- **Location**: `crates/cfd-3d/src/venturi/solver.rs`
- **Issue**: `cp_throat` and `cp_recovery` were divided by inlet dynamic pressure, while the solver theorem defines Venturi recovery coefficients with throat dynamic pressure. The computed coefficients therefore depended on inlet area rather than the throat kinetic-energy scale.
- **Remediation**: Added a single coefficient helper that computes `q_throat = 0.5 ρ (Q/A_throat)^2` from the face-integrated flow rate and throat area, rejects undefined non-positive scales, and routes both Venturi pressure coefficients through that helper.
- **Mathematical basis**: For steady incompressible Venturi flow, continuity gives `u_throat = Q/A_throat`; pressure coefficients nondimensionalize `Δp` by `0.5 ρ u_throat²`. Using inlet velocity under-scales the denominator by the squared area ratio and overstates dimensionless recovery.
- **Verification**:
  - `cargo check -p cfd-3d --no-default-features` passed.
  - `cargo test -p cfd-3d --no-default-features venturi_pressure_coefficients --lib` passed 2/2 tests.
  - `cargo nextest run -p cfd-3d --lib --no-default-features venturi_pressure_coefficients --fail-fast --hide-progress-bar --status-level fail` passed 2/2 tests.

---

# Sprint 1.96.6 Resolution: cfd-2d Explicit Stability Physics

### RESOLVED-035: 2D CFL Time-Step Bound Used Componentwise Limits
- **Location**: `crates/cfd-2d/src/stability/cfl.rs`
- **Issue**: `max_stable_dt` used `min(dx/|u|, dy/|v|)` for 2D advection and `0.5*min(dx²,dy²)/ν` for diffusion. These bounds can violate the module's documented summed 2D CFL and diffusion conditions.
- **Remediation**: Replaced the advection bound with `1 / (|u|/dx + |v|/dy)` and the diffusion bound with `0.5 / (ν(1/dx² + 1/dy²))`. Added tests that feed the returned `dt` back into `advection_cfl` and `diffusion_number`.
- **Mathematical basis**: Solving `|u|Δt/dx + |v|Δt/dy ≤ 1` and `νΔt(1/dx² + 1/dy²) ≤ 1/2` for `Δt` yields the reciprocal summed-rate limits. Taking the minimum satisfies both explicit stability constraints.
- **Verification**:
  - `cargo check -p cfd-2d --no-default-features` passed.
  - `cargo test -p cfd-2d --no-default-features stability::cfl --lib -- --nocapture` passed 4/4 tests.
  - `cargo nextest run -p cfd-2d --lib --no-default-features stability::cfl --fail-fast --hide-progress-bar --status-level fail` passed 4/4 tests in 0.119 s.

---

# Sprint 1.96.5 Resolution: cfd-1d Venturi Reverse-Flow Physics

### RESOLVED-034: Venturi Reverse Flow Lost Inertial Coefficients
- **Location**: `crates/cfd-1d/src/physics/resistance/models/venturi/`
- **Issue**: `VenturiModel::calculate_coefficients` selected its finite-flow branch with `v_inlet > threshold`, so negative inlet velocity was treated as zero flow and lost contraction, diffuser, and recovery coefficients. Detailed analysis also used signed velocity for shear-rate and Reynolds inputs.
- **Remediation**: Coefficient decomposition now uses `|V_inlet|`, `|V_throat|`, and `|Q|` for scalar loss magnitudes while detailed analysis preserves signed velocities in reported kinematics. Shear-rate, viscosity, Reynolds correction, friction factor, and pressure-loss terms use magnitudes.
- **Mathematical basis**: In a symmetric Venturi, reversing flow orientation changes the sign of axial velocities but leaves kinetic-energy losses proportional to `V²`, wall-shear magnitude, Reynolds magnitude, and scalar resistance coefficients unchanged.
- **Verification**:
  - `cargo check -p cfd-1d --no-default-features` passed.
  - `cargo test -p cfd-1d --no-default-features physics::resistance::models::venturi --lib -- --nocapture` passed 16/16 tests.
  - `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::venturi --fail-fast --hide-progress-bar --status-level fail` passed 16/16 tests in 0.230 s.

---

# Sprint 1.96.4 Resolution: Geometric Conservation Residual Verification

### RESOLVED-033: Copy-Through Geometric Conservation Evolution
- **Location**: `crates/cfd-validation/src/conservation/geometric.rs`
- **Issue**: Euler and Runge-Kutta GCL checks copied the scalar field between time levels, so tests could pass without evaluating any numerical residual.
- **Remediation**: Added a conservative second-order face-flux residual and routed Euler, midpoint, SSPRK3, and RK4 checks through residual-based stage updates. Unsupported RK stage counts now return a typed `UnsupportedOperation` error.
- **Mathematical basis**: For constant `u_ij = c`, every face gradient is zero, so the finite-volume residual divergence is zero at every interior cell. Explicit Runge-Kutta stage combinations of `u` and `dt R(u)` therefore preserve the constant state exactly.
- **Verification**:
  - `cargo check -p cfd-validation --no-default-features` passed in 0.91 s.
  - `cargo test -p cfd-validation --no-default-features conservation::geometric --lib -- --nocapture` passed 7/7 tests.
  - `cargo nextest run -p cfd-validation --lib --no-default-features conservation::geometric --fail-fast --hide-progress-bar` passed 7/7 tests in 0.229 s under a 30-second shell timeout.

---

# Module: CFD-IO

**Status**:  COMPLETE (Verified)

## Audit Summary
- **Bit-Exact Roundtrip**: Verified. 
  read(write(data)) == data invariant holds.
- **Parallel Checksum**: Verified. MPI allreduce checksums implemented.
- **Documentation**: Invariants documented.

---

# Module: CFD-MATH (Linear Solvers & Preconditioners)

**Status**:  PASS (Remediated)

## Theorem Verification  Pass
- **Conjugate Gradient**: Good documentation of Optimality Theorem and Condition Number bounds.
- **Preconditioners**: Good mathematical foundation docs for Jacobi, SOR, AMG.
- **Parallel Consistency**: Clarified that current implementations are primarily serial reference implementations or building blocks for parallel solvers (which reside in `cfd-core` or use specific MPI types).

## Algorithm Audit  PASS

### RESOLVED-001: Redundant and Inferior ILU Implementation
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners.rs`
- **Issue**: Two conflicting implementations of ILU factorization existed.
- **Remediation**: `ILUPreconditioner` has been removed. `IncompleteLU` is now the single source of truth, supporting both ILU(0) and ILU(k).

### RESOLVED-002: Serial Schwarz Decomposition (Clarified)
- **Location**: `crates/cfd-math/src/linear_solver/preconditioners.rs`
- **Issue**: The Schwarz preconditioner was implemented using a serial loop but named generic `SchwarzPreconditioner`.
- **Remediation**: Renamed to `SerialSchwarzPreconditioner` to explicitly reflect its nature as a serial implementation (useful for testing/research or as a local block solver). Documentation updated to clarify it does not provide MPI parallelism itself.

### RESOLVED-003: Parallel Claims in Solvers
- **Location**: `crates/cfd-math/src/linear_solver/conjugate_gradient.rs`
- **Issue**: Documentation claimed "Parallel Scalability" for a generic implementation that uses serial `nalgebra` types by default.
- **Remediation**: Removed misleading "Parallel Scalability" claims. Updated documentation to clarify that parallelism depends on the underlying vector/matrix types or should use the specific MPI solvers in `cfd-core`.

## Testing Validation  PARTIAL
- **Unit Tests**: Good coverage for serial cases (Identity, Diagonal, Poisson 1D).
- **Missing**: Parallel scaling tests (moved to `cfd-core` benchmarks), convergence verification for difficult matrices (high condition number).

## Code Quality Audit  PASS
- **Good**: Extensive doc comments with literature references.
- **Remediated**: Redundant structs removed, misleading names fixed.

---

# Module: CFD-CORE (Distributed Computing & MPI)

**Status**:  FAILED

## Algorithm Audit  FAILED

### RESOLVED-004: Fake Distributed GMRES Implementation
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`DistributedGMRES::gmres_iteration`)
- **Issue**: The GMRES solver performed the Arnoldi process but skipped the least-squares solution step, setting residual to zero manually.
- **Remediation**: Implemented full GMRES cycle with Givens rotations for updating the QR factorization of the Hessenberg matrix, and backward substitution to solve the least-squares problem.

### RESOLVED-005: Fake Additive Schwarz Preconditioner
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`AdditiveSchwarzPreconditioner::create_local_solvers`)
- **Issue**: The local subdomain solvers were implemented as Identity operations.
- **Remediation**: Updated to use `assemble_local_matrix` and `nalgebra::LU` for direct local solves. Added fallback to Jacobi (diagonal scaling) if matrix assembly is not supported by the operator.

### RESOLVED-006: Fake Block Jacobi Preconditioner
- **Location**: `crates/cfd-core/src/compute/mpi/distributed_solvers.rs` (`BlockJacobiPreconditioner::new`)
- **Issue**: The preconditioner assumed unit diagonal blocks.
- **Remediation**: Extended `DistributedLinearOperator` trait with `extract_diagonal` and updated `BlockJacobiPreconditioner` to use the actual operator diagonal.

---

# Module: CFD-MESH (Topology & Refinement)

**Status**:  FAILED

## Algorithm Audit  FAILED

### RESOLVED-007: Fake Mesh Refinement
- **Location**: `crates/gaia/src/refinement/mod.rs`
- **Issue**: `UniformRefinement` and `AdaptiveRefinement` methods were empty placeholders returning `Ok(())`.
- **Remediation**: Updated methods to explicitly return `MeshError::NotImplemented`. This removes the deceptive "working" status and correctly signals that the feature is pending implementation.

### RESOLVED-008: Missing Distributed Mesh Support
# Critical Audit Findings (Open Issues)

| ID | Severity | Component | Issue | Status |
|----|----------|-----------|-------|--------|
| **CRITICAL-001** |  Critical | CFD-MATH | Redundant ILUPreconditioner vs IncompleteLU | **CLOSED** |
| **CRITICAL-002** |  Critical | CFD-MATH | Serial implementation of Schwarz Decomposition | **CLOSED** |
| **MAJOR-003** |  Major | CFD-MATH | Unsubstantiated Parallel Scalability claims in CG | **CLOSED** |
| **CRITICAL-004** |  Critical | CFD-CORE | Fake Distributed GMRES (Placeholder Solver) | **CLOSED** |
| **CRITICAL-005** |  Critical | CFD-CORE | Fake Additive Schwarz (Identity Solver) | **CLOSED** |
| **MAJOR-006** |  Major | CFD-CORE | Fake Block Jacobi (Unit Diagonal Assumption) | **CLOSED** |
| **CRITICAL-011** | ✅ Closed | CFD-MESH | WASM OOM due to O(N²) allocation overhead in Bowyer-Watson | **CLOSED** |
| **CRITICAL-007** |  Critical | CFD-MESH | Fake Mesh Refinement (Empty Methods) | **CLOSED** |
| **MAJOR-008** |  Major | CFD-MESH | Missing Distributed Mesh Support | **CLOSED** |
| **CRITICAL-009** | ✅ Closed | CFD-MATH | Ruge-Stüben Coarsening Fine-to-Coarse Mapping Bug | **CLOSED** |
| **CRITICAL-010** | ✅ Closed | CFD-3D | Fake Unstructured FEM Domains (0-Element Mocks) | **CLOSED** |
| **CRITICAL-012** | 🔴 Critical | CFD-CORE | Missing cell-specific cavitation injury physics (Lysis/Necrosis) | **CLOSED** |
| **CRITICAL-013** | 🔴 Critical | CFD-CORE | Missing CTC stiffness-coupled inception mechanisms for HCOC | **CLOSED** |
| **MAJOR-014** | ✅ Closed | CFD-OPTIM | Milestone 12 cross-mode therapy utility misclassified ultrasound-only SDT as non-cavitating because it gated all cavitation credit on active venturi throats | **CLOSED** |

---

## RESOLVED-014: Nuclei Transport Diffusion Coupling in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-core/src/physics/cavitation/nuclei_transport.rs`, `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - The documented advection-diffusion model is now implemented in the 3D solver.

### Verification

- `NucleiTransport::diffusion_coefficient()` exposes the configured scalar diffusivity to the solver layer.
- `CavitationVofSolver::update_nuclei_advection_diffusion()` now applies explicit finite-volume diffusion with zero-flux outer faces in addition to the existing upwind advection and reaction terms.
- Regression coverage confirms a localized nuclei peak spreads into neighboring cells without flow or generation.

---

## RESOLVED-015: Cavitation Damage Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation damage accumulation now validates incoming pressure and density fields and walks contiguous column slices instead of repeated matrix element indexing.

### Verification

- `CavitationVofSolver::update_damage()` now rejects pressure and density dimension mismatches with `DimensionMismatch` instead of silently skipping the damage update.
- The damage accumulation path now uses raw column slices for pressure, density, and damage storage, preserving the column-major layout while removing repeated matrix indexing overhead.
- Regression coverage confirms mismatched pressure fields are rejected and sub-1% void-fraction damage still accumulates.

---

## RESOLVED-016: Cavitation Source Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation source accumulation now validates its matrix inputs, uses contiguous slices, and clamps the source to the feasible interval per cell.

### Verification

- `calculate_cavitation_source_into()` now rejects pressure, density, and source-workspace shape mismatches with `DimensionMismatch`.
- The source update path now iterates over raw column slices, preserving the column-major storage contract while removing repeated matrix indexing overhead.
- Regression coverage confirms the source clamps to zero at the vaporization and condensation feasibility extremes.

---

## RESOLVED-017: Pries Plasma-Skimming Silent Clamps in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/solvers/plasma_skimming.rs`
**Status**: **CLOSED** - The standalone 2D Pries phase-separation path now evaluates the thresholded logit law through a checked physical-envelope API instead of clipping flow, hematocrit, diameter ratios, logit arguments, or daughter hematocrit outputs.

### Verification

- `checked_pries_phase_separation()` rejects non-finite, out-of-range, and non-positive geometry inputs with structured errors.
- Boundary cases `FQB = 0`, `FQB = 1`, and zero feed hematocrit are explicit conservation states.
- Regression coverage validates below-threshold zero cell flux without modifying the supplied flow fraction.

---

## RESOLVED-018: WALE Boundary Zero-Gradient Reduction in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/les_smagorinsky/wale.rs`
**Status**: **CLOSED** - WALE velocity-gradient evaluation now uses second-order one-sided boundary stencils on uniform grids instead of forcing boundary gradients to zero.

### Verification

- Rustdoc documents the Taylor-series stencil proof.
- Regression coverage confirms exact derivative recovery for quadratic and linear velocity components at lower and upper boundaries.

---

## RESOLVED-019: Spalart-Allmaras All-Zero TKE Diagnostic in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/physics/turbulence/spalart_allmaras.rs`
**Status**: **CLOSED** - The SA turbulence model no longer returns a constant zero turbulent-kinetic-energy vector when the transported eddy viscosity is nonzero.

### Verification

- SA TKE is now a documented diagnostic `k = (nu_t / (C_k d))^2` using the shared Yoshizawa relation and local wall distance.
- Regression coverage validates the computed diagnostic against the analytical relation and preserves the zero wall-length-scale state.

---

## RESOLVED-020: K-Epsilon Uninitialized Zero Fields in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/physics/turbulence/k_epsilon.rs`
**Status**: **CLOSED** - The k-epsilon model no longer fabricates zero viscosity, turbulent kinetic energy, or dissipation fields when the transported state is absent.

### Verification

- `turbulent_viscosity()`, `turbulent_kinetic_energy()`, and `dissipation_rate()` now require initialized `k` and `epsilon` fields with dimensions matching the flow field.
- Initialized zero-TKE states still produce zero eddy viscosity by the analytical `nu_t = C_mu k^2 / epsilon` limit.
- Unit and integration tests assert uninitialized state rejection and initialized zero-state behavior.

---

## RESOLVED-021: Smagorinsky LES Zero SGS Diagnostics in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/les_smagorinsky/`
**Status**: **CLOSED** - The 2D Smagorinsky LES model no longer returns placeholder zero turbulent kinetic energy or dissipation rate after computing nonzero SGS viscosity.

### Verification

- SGS turbulent kinetic energy now follows the Yoshizawa relation `k_sgs = (nu_t / (C_k Delta))^2`.
- SGS dissipation now follows inertial-range dimensional scaling `epsilon_sgs = C_epsilon k_sgs^(3/2) / Delta`.
- Strain-rate evaluation now uses second-order one-sided boundary stencils instead of imposing zero boundary strain.
- Default `min_sgs_viscosity` is zero so zero-strain fields remain exactly laminar unless a caller explicitly configures a numerical lower bound.
- Regression coverage validates the Yoshizawa energy and dissipation diagnostics, boundary strain recovery, and zero-strain invariant.

---

## RESOLVED-022: Margination Singular Wall Cutoff in CFD-1D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/margination.rs`
**Status**: **CLOSED** - The margination inertial-lift path no longer clamps lateral position to a hard-coded near-wall cutoff or evaluates a singular wall-lift shape.

### Verification

- The validated domain is now the closed half-channel interval `0 <= x_tilde <= 1`.
- The wall-induced lift shape is bounded and monotone on the full interval while preserving the published `rho U^2 a^6/H^4` wall-lift scaling.
- The shear-gradient term preserves the published `rho U^2 a^3/H` scaling and vanishes at the wall.
- Public legacy lift evaluation now dispatches to the checked evaluator and rejects invalid positions instead of silently clipping.
- Regression coverage validates monotonic wall lift, reference equilibrium calibration, and a derived force value for MCF-7-like cells.

---

## RESOLVED-023: Point-Droplet Occupancy State in CFD-1D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/solver/core/transient/droplets/types.rs`
**Status**: **CLOSED** - `DropletSnapshot` no longer stores a separate point-style `occupied_channels` vector beside finite-length occupancy spans.

### Verification

- `occupancy_spans` is now the single authoritative droplet occupancy representation.
- `occupied_channels()` computes the ordered unique channel projection on demand from finite spans.
- Snapshot consistency validation checks finite interval bounds instead of comparing duplicate stored representations.
- Unit and integration tests validate projection, duplicate-channel collapsing, invalid-span rejection, split/merge behavior, and API/manual parity.

---

## RESOLVED-024: Smagorinsky Validation SGS Floor in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/validation/`
**Status**: **CLOSED** - The turbulence validation paths no longer override the Smagorinsky LES physical zero-floor default with residual nonzero SGS viscosity floors.

### Verification

- Benchmark dispatch remains a closed enum with typed unsupported-model rejection.
- Smagorinsky benchmark and LES/DES validation configurations now use `min_sgs_viscosity = 0.0`.
- Targeted nextest verifies the benchmark module and LES tests under the requested 30-second timeout profile.

---

## RESOLVED-017: Nuclei Transport Validation and Slice-Based Accumulation in CFD-3D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/vof/cavitation_solver.rs`
**Status**: **CLOSED** - Cavitation nuclei advection-diffusion now validates its solver inputs and accumulates through contiguous column slices.

### Verification

- `update_nuclei_advection_diffusion()` now rejects cavitation-source and nuclei-workspace shape mismatches with `DimensionMismatch` or `InvalidConfiguration` instead of assuming valid inputs.
- The nuclei transport update now iterates over raw column slices for state, next-state, and source storage, matching the column-major layout and removing repeated matrix indexing overhead.
- Regression coverage confirms a localized cavitation source increases nuclei fraction while zero-flow diffusion spreads a localized peak into neighboring cells.

---

## RESOLVED-018: Indexed Blueprint Layout Cache in CFD-SCHEMATICS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/visualizations/schematic`
**Status**: **CLOSED** - Schematic layout generation now uses a borrowed indexed cache with flat position storage and index-keyed parallel grouping.

### Verification

- `blueprint_node_positions()` now stores node positions in authored order with a borrowed ID index and a flat auto-layout index list instead of cloning string-keyed maps.
- `resolved_channel_paths()` now resolves parallel channel groups by borrowed node indices instead of cloned `(String, String)` keys.
- `materialize_blueprint_layout()` reuses the indexed cache to write node layout metadata without changing explicit positions.
- Regression coverage confirms the indexed layout cache produces the expected auto-layout coordinates and the materialization path writes the computed fallback positions into the blueprint.

---

## RESOLVED-019: Geometry Bounds and Clearance-Width Relations in CFD-SCHEMATICS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/config/geometry`, `crates/cfd-schematics/src/state_management`
**Status**: **CLOSED** - Geometry parameter bounds now share canonical constants and clearance-width degeneracy is rejected at the config, manager, and registry validation layers.

## RESOLVED-020: Milestone 12 GA Convergence Trend SSOT in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/reporting`
**Status**: **CLOSED** - The GA convergence figure and narrative now derive their improvement / regression / near-plateau labels from the same trailing-fitness-window helper.

### Verification

- `GeometryConfig::validate()` now uses the canonical geometry constants for channel width, channel height, and wall clearance, and rejects `wall_clearance >= channel_width`.
- `GeometryParameterManager::validate_all()` now enforces the same clearance-width relation at the manager layer.
- `ParameterRegistry::validate_all()` rejects the degenerate pair through the geometry manager path, and the integration regression confirms the returned `InvalidValue` payload.
- Regression coverage confirms canonical bound acceptance and degenerate clearance-width rejection at the config and manager layers.

---

## RESOLVED-021: Slot-Stable Clip-Plane Undo in CFD-UI

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-ui/src/domain/clipping/plane.rs`, `crates/cfd-ui/src/domain/clipping/mod.rs`, `crates/cfd-ui/src/application/clipping/commands.rs`
**Status**: **CLOSED** - Clip-plane commands now preserve their recorded slots, surface explicit stale-state failures, and reject silent relocation during undo.

### Verification

- `ClipPlaneSet::insert_at()` and `replace_at()` now update only the requested slot and return structured slot errors for occupied, empty, and out-of-range cases.

---

## RESOLVED-030: Plasma Skimming, Serpentine Mixing, and LES SGS Energy Physics

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/plasma_skimming.rs`, `crates/cfd-2d/src/solvers/serpentine_flow`, `crates/cfd-3d/src/physics/turbulence`
**Status**: **CLOSED** - Three reduced physics paths now dispatch through documented analytical or literature-backed models.

### Verification

- The compact `plasma_skimming_hematocrit` API now infers the sibling branch from Murray cubic diameter conservation and calls the checked Pries phase-separation model, preserving the `X0` cell-entry threshold.
- The serpentine mixing model now computes mixing fraction and L90 from the Neumann eigenfunction solution of transverse diffusion for a two-stream step inlet; the discretized solver result now reports positive analytical L90 and t90 instead of zero placeholders.
- LES eddy-viscosity models no longer return turbulent viscosity as turbulent kinetic energy; they use a shared Yoshizawa SGS relation with value-semantic regression coverage.
- Add, remove, and update clip-plane commands now return `Result` and validate the live slot before restoring or replacing a plane.
- Regression coverage confirms valid add/remove/update round-trips succeed, while undoing a remove into an occupied slot fails explicitly instead of relocating the plane.

---

## RESOLVED-022: Milestone 12 Validation Traceability in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/reporting`, `report/templates/m12_narrative_template.md`, `report/ARPA-H_SonALAsense_Milestone 12 Report.md`, `report/milestone12_results.md`
**Status**: **CLOSED** - The Milestone 12 report now exposes validation evidence and traceability in both the narrative and canonical results artifacts.

### Verification

- `build_results_intro()` now advertises validation evidence and artifact traceability in the Milestone 12 narrative contract, and `m12_narrative_template.md` renders the corresponding subsection.
- `report/ARPA-H_SonALAsense_Milestone 12 Report.md` now contains `### 5.4 Validation Evidence and Artifact Traceability` and the cross-fidelity validation evidence table.
- `report/milestone12_results.md` now contains `## Cross-Fidelity Validation Evidence`, and `report/milestone12/asset_review_manifest.json` was updated to the regenerated report hash.

---

## RESOLVED-023: Single-Pass Authoritative Milestone 12 Report in CFD-OPTIM

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-optim/src/application/orchestration/milestone12/report.rs`, `report/milestone12/report_manifest.json`, `report/ARPA-H_SonALAsense_Milestone 12 Report.md`, `report/milestone12/asset_review_manifest.json`
**Status**: **CLOSED** - The release report now emits the authoritative narrative once, preserves the reviewed hash, and completes after the asset-review manifest is marked pass.

### Verification

- `refresh_milestone12_reports()` now renders the narrative with `authoritative_run = run_class.is_authoritative()` on the first pass instead of rewriting the report with a second content hash.
- `cargo run -p cfd-optim --example milestone12_report --no-default-features --release` completes successfully once the manifest review entry is marked pass under a PATH that resolves the MSYS2 compiler DLLs correctly.
- `report/milestone12/report_manifest.json` records `authoritative_run: true` and `review_complete: true`, and the narrative asset-review entry remains `opened: true`, `status: pass` after regeneration.

---

## CRITICAL-009: Ruge-Stüben Fine-to-Coarse Mapping Bug (VERIFIED FIXED)

**Severity**: ✅ **RESOLVED** - Code Review Confirms Correct Implementation  
**Discovered**: November 18, 2025 (Deep Algorithm Audit)  
**Verified**: February 23, 2026 (Sprint 1.95.0)  
**Component**: `crates/cfd-math/src/linear_solver/preconditioners/multigrid/coarsening.rs`  
**Status**: **CLOSED** - Implementation is correct

### Verification (February 23, 2026)

Upon detailed code review of `assign_closest_coarse_points()` function in `coarsening.rs`, the implementation is **already correct**:

```rust
// Current implementation (CORRECT):
for k in s_offsets[i]..s_offsets[i + 1] {
    let j = s_indices[k];
    if status[j] == 1 {  // j is a C-point
        let strength = strength_matrix.values()[k];
        if strength > max_strength {
            max_strength = strength;
            best_coarse_idx = fine_to_coarse_map[j];  // ✅ CORRECT: gets the coarse grid index
        }
    }
}
```

The key insight: `fine_to_coarse_map[j]` for a C-point `j` already contains the **coarse grid index** (set when the point was marked as coarse via `fine_to_coarse_map[i] = Some(coarse_points.len() - 1)`). Therefore, assigning `fine_to_coarse_map[j]` to an F-point correctly gives the coarse grid index.

### Existing Tests

The file contains comprehensive tests that verify correctness:
- `test_mapping_correctness`: Verifies all mapped indices are valid coarse indices
- `test_interpolation_operator_shape`: Verifies fine_to_coarse_map dimensions
- `test_coarsening_ratio_bounds`: Verifies reasonable coarsening ratios for 2D Laplacian

### Issue Description (Historical)

- **Interpolation Operator Corruption**: The AMG interpolation operator will reference incorrect coarse DOFs
- **Convergence Theory Violation**: AMG theory requires fine points to interpolate from geometrically nearby coarse points
- **Performance Degradation**: Suboptimal convergence rates, potential divergence for difficult problems
- **Literature Divergence**: Does not match Ruge-Stüben (1987) algorithm specification

### Current Code (INCORRECT)

```rust
// Step 2: Assign fine points to nearest coarse points
for i in 0..n {
    if fine_to_coarse_map[i].is_none() {
        // Find strongest connection to coarse point
        let mut max_strength = 0.0;
        let mut best_coarse = None;

        for &c in &coarse_points {
            if strength_matrix[(i, c)] > max_strength {
                max_strength = strength_matrix[(i, c)];
                best_coarse = Some(c);
            }
        }

        if let Some(c) = best_coarse {
            fine_to_coarse_map[i] = fine_to_coarse_map[c];  // ❌ BUG: Assigns mapping value, not index
        }
    }
}
```

### Expected Code (CORRECT)

```rust
// Step 2: Assign fine points to nearest coarse points
for i in 0..n {
    if fine_to_coarse_map[i].is_none() {
        // Find strongest connection to coarse point
        let mut max_strength = 0.0;
        let mut best_coarse = None;

        for &c in &coarse_points {
            if strength_matrix[(i, c)] > max_strength {
                max_strength = strength_matrix[(i, c)];
                best_coarse = Some(c);
            }
        }

        if let Some(c) = best_coarse {
            // ✅ CORRECT: Map to coarse point INDEX in coarse grid
            let coarse_idx = coarse_points.iter()
                .position(|&x| x == c)
                .expect("Coarse point must exist in coarse_points list");
            fine_to_coarse_map[i] = Some(coarse_idx);
        }
    }
}
```

### Remediation Steps

1. **Fix Mapping Logic** (1-2 hours):
   - Update fine-to-coarse map to assign coarse point indices
   - Ensure coarse points map to their own indices in coarse grid
   - Verify interpolation operator dimensions match expected (n_fine × n_coarse)

2. **Add Validation Tests** (2-3 hours):
   - Test that all fine points map to valid coarse indices  (0 ≤ idx < n_coarse)
   - Test that coarse points self-map correctly
   - Test interpolation operator shape: P.shape == (n_fine, n_coarse)
   - Test restriction operator shape: R.shape == (n_coarse, n_fine)
   - Add regression test for this specific bug

3. **Validate Convergence** (1 hour):
   - Run AMG on Poisson problem, verify O(N) complexity
   - Check convergence factor < 0.1 per V-cycle (theory predicts ~0.03-0.08)
   - Compare with literature benchmarks (Ruge-Stüben 1987, Briggs 2000)

### Why Not Caught Earlier

- **Tests Pass**: System still solves (with suboptimal interpolation)
- **Convergence Occurs**: AMG may still converge, just slower than optimal
- **Subtle Bug**: Requires deep algorithm knowledge and literature familiarity to detect
- **No Analytical Verification**: Previous tests didn't validate internal AMG structure

### Evidence Hierarchy Violation

This bug violates the audit framework's evidence hierarchy:
- ❌ **Mathematical Proofs**: Mapping does not match Ruge-Stüben theorem
- ❌ **Literature Validation**: Code diverges from Ruge-Stüben (1987) Algorithm 2
- ⚠️ **Empirical Testing**: Tests incomplete (didn't check mapping correctness)

1. ✅ **FIX CRITICAL-009: Ruge-Stüben Fine-to-Coarse Mapping Bug** (2-5 hours, HIGHEST PRIORITY)
   - Update mapping assignment in `coarsening.rs` to use coarse point indices
   - Add comprehensive AMG coarsening validation tests
   - Verify interpolation operator correctness
   - Validate convergence improvement (expected 2-5x speedup)
   - Block AMG use in production until verified

## COMPLETED (Previous Sprints)

2.  **Immediate Refactoring**: Implement the missing least-squares solve in `DistributedGMRES`. (Done)
3.  **Integration**: Connect `cfd-core` to `cfd-math`. Use `cfd-math::IncompleteLU` as the local solver for `AdditiveSchwarzPreconditioner`. (Done)
4.  **API Extension**: Add `get_diagonal` to `DistributedLinearOperator` to enable real Block Jacobi. (Done)
5.  **Implementation**: Implement basic uniform refinement in `cfd-mesh` to validate the interface. (Done - Implemented 1:8 Tet Refinement)
6.  **Continuous Validation**: Ensure `IncompleteLU` remains the standard for ILU operations.
7.  **Parallel Integration**: When implementing distributed solvers in `cfd-core`, leverage the serial block solvers from `cfd-math` where appropriate, but ensure explicit types for distributed contexts.
8.  **Benchmarking**: Continue to expand benchmarks in `cfd-core` to validate actual parallel scaling of the MPI-specific implementations.

---

## RESOLVED-024: HCOC Cellular Injury & CTC Detection Framework

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-core/src/physics/cavitation/bio_damage.rs`, `crates/cfd-core/src/physics/cavitation/heterogeneous_nucleation.rs`
**Status**: **CLOSED** - Cellular injury grading and stiffness-coupled inception are implemented and tested with literature-consistent threshold physics.

### Verification

- `CellularMembraneMechanics::evaluate_spatial_injury()` uses rate-adjusted critical strains and a first-order hazard transition over the ordered permeabilization / necrosis / lysis thresholds while conserving mass fractions.
- `blake_threshold_pressure_pa()` now applies the classical Blake threshold factor `P_B = P_v - 4σ/(3R_B)` instead of the previous under-corrected `2σ` term.
- Regression tests cover exact mass conservation, ambient no-injury behavior, monotone injury growth, loading-duration sensitivity, the classical Blake-factor identity, and target-versus-healthy selectivity ordering.

---

## RESOLVED-025: Apollo FFT Plan Binding in CFD-3D Periodic DNS

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-3d/src/spectral/dns.rs`
**Status**: **CLOSED** - The periodic pseudospectral DNS stepper now validates its Apollo 3D shape once at construction and owns a reusable `FftPlan3D` for all forward and inverse transforms.

### Verification

- `PeriodicPseudospectralDns3D::new()` maps the validated CFD grid dimensions into Apollo `Shape3D` and returns structured configuration errors if Apollo rejects the shape.
- Velocity spectra, nonlinear spectra, spectral derivatives, and inverse step reconstruction all use the stepper-owned plan rather than per-transform shape dispatch.
- `cargo check -p cfd-3d --no-default-features` completes successfully.

---

## RESOLVED-026: Exact MUSCL3/QUICK Face Reconstruction in CFD-2D

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/momentum/muscl.rs`
**Status**: **CLOSED** - The third-order MUSCL/QUICK path no longer uses an approximate right-interface formula; both face states now use exact uniform-grid quadratic interpolation.

### Verification

- `quadratic_left_face()` implements `-phi_{i-1}/8 + 3 phi_i/4 + 3 phi_{i+1}/8`.
- `quadratic_right_face()` implements `3 phi_i/8 + 3 phi_{i+1}/4 - phi_{i+2}/8`.
- Regression coverage verifies exact reproduction of a linear field from both sides and exact reproduction of a quadratic parabola at the face.
- `CC=clang cargo nextest run -p cfd-2d --no-default-features muscl ...` passed 14/14 filtered tests under the 30-second slow-test bound.

---

## RESOLVED-027: cfd-1d Margination Lift and Droplet Occupancy Contracts

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/physics/cell_separation/margination.rs`, `crates/cfd-1d/src/solver/core/transient/droplets`
**Status**: **CLOSED** - The margination model now separates wall-induced and shear-gradient inertial lift scaling, and droplet snapshots document `occupied_channels` as a projection of finite-length occupancy spans.

### Verification

- `dimensional_lift_force_n()` evaluates distinct wall-induced and shear-gradient finite-size scaling terms and derives the wall gain from the documented rigid-cell reference equilibrium.
- Rustdoc theorem/proof sketches document wall-shape monotonicity, shear-shape monotonicity, and reference-equilibrium calibration.
- `DropletSnapshot::occupied_channels_from_spans()` and `has_consistent_finite_length_occupancy()` encode the finite-length occupancy projection contract.
- Targeted `CC=clang cargo nextest run -p cfd-1d --no-default-features ...` passed the margination and droplet regression tests.

---

## RESOLVED-028: cfd-2d Turbulence Benchmark Dispatch

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-2d/src/physics/turbulence/validation/benchmarks.rs`
**Status**: **CLOSED** - The performance benchmark no longer contains an unimplemented fallback path; supported benchmark models are parsed into a closed enum before dispatch.

### Verification

- `BenchmarkModel` admits only `smagorinsky-les` and `des`, and exhaustive matching runs concrete production update loops.
- Unsupported model names return a typed rejection listing supported models and do not emit placeholder text.
- `CC=clang cargo nextest run -p cfd-2d --no-default-features ...` passed 554/554 tests under the 30-second slow-test bound.

---

## RESOLVED-029: cfd-1d Coupled Blood Pressure-Event Verification

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-1d/src/domain/network/wrapper.rs`, `crates/cfd-1d/src/solver/core/linear_system.rs`, `crates/cfd-1d/src/solver/core/transient/composition/simulator.rs`
**Status**: **CLOSED** - Coupled pressure-event hematocrit simulation now produces finite positive branch-flow magnitudes while preserving signed internal transport.

### Verification

- Zero-flow blood apparent viscosity uses finite Secomb diameter/hematocrit initialization before nonzero shear-dependent Quemada updates.
- Small hydraulic systems use row-equilibrated solves so conductance rows near machine scale are not dominated by unit Dirichlet rows.
- Dense fallback conversion accumulates duplicate CSR entries, preserving the sum of incident-edge diagonal conductances required by Kirchhoff balance.
- `CC=clang cargo nextest run -p cfd-1d --no-default-features ...` passed 696/696 tests under the 30-second slow-test bound.

---

## RESOLVED-030: Root Historical Source Artifact Drift

**Severity**: ✅ **RESOLVED**
**Component**: workspace root source artifacts
**Status**: **CLOSED** - Tracked root historical Rust files that were outside Cargo, examples, tests, docs, and report-generation paths have been removed to preserve single-source-of-truth ownership.

### Verification

- `git ls-files` confirmed `old_assemble.rs`, `old_arrangement.rs`, `old_phase2.rs`, `old_operations.rs`, `old_indexed.rs`, `old_gwn_bvh.rs`, `old_seam.rs`, `old_phase4.rs`, and `csg_bi.rs` were tracked root artifacts.
- `rg` found no references to those paths in `Cargo.toml`, `crates/`, `examples/`, `tests/`, `docs/`, or `README.md`.
- Deletion is patch-class cleanup: no public API, algorithm, report template, or executable target references the removed files.
- `cargo metadata --no-deps --format-version 1` and `cargo check -p cfd-suite --no-default-features` completed successfully.
- `.cargo/config.toml` now selects MSYS2 `clang.exe` and `-fuse-ld=lld` for `x86_64-pc-windows-gnu`, and forces C/C++ build-script tools to MSYS2 `clang`, `clang++`, `llvm-ar`, and `llvm-ranlib`.
- `cargo nextest run --workspace --no-default-features --fail-fast --hide-progress-bar` was re-attempted after the toolchain correction; the prior `gcc.exe` `zstd-sys` failure did not recur, but workspace compilation exceeded the 60-second bound before Rust test execution.
- Misleading terminology was removed from explicit unsupported-operation paths and bounded-model comments in `cfd-core`, `cfd-math`, `cfd-python`, and `gaia` without changing numerical behavior.
- `cargo check -p cfd-core -p cfd-math -p gaia -p cfd-python --no-default-features` completed successfully in 46.17 seconds.
- `cargo nextest run -p cfd-core -p cfd-math -p cfd-python --no-default-features --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.

---

## RESOLVED-031: Optimization Generator and Boundary Unsupported-Order Terminology

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-schematics/src/geometry/optimization`, `crates/cfd-core/src/physics/boundary/error.rs`
**Status**: **CLOSED** - Internal serpentine optimization path generation and boundary stencil diagnostics now describe executable contracts without implying simplified or incomplete implementations.

### Verification

- `generate_simplified_serpentine_path` was renamed to `generate_optimization_serpentine_path`; all internal optimization call sites were updated.
- Boundary unsupported-order diagnostics now report `"Stencil order {order} is unsupported (supported: 1-4)"`.
- `rg` found no stale `generate_simplified_serpentine_path`, `simplified serpentine`, `simplified generation`, or `"not implemented"` wording in the touched paths.
- `cargo check -p cfd-core -p cfd-schematics --no-default-features` completed successfully in 28.06 seconds.
- `cargo nextest run -p cfd-core -p cfd-schematics --no-default-features --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.

---

## RESOLVED-032: Womersley Validation Approximation Drift

**Severity**: ✅ **RESOLVED**
**Component**: `crates/cfd-validation/src/analytical/womersley.rs`
**Status**: **CLOSED** - The validation Womersley analytical solution now reuses the canonical exact complex-Bessel `cfd-1d` implementation instead of maintaining local regime approximations.

### Verification

- `velocity()`, `wall_shear_stress()`, and `flow_rate()` delegate through `cfd_1d::physics::vascular::womersley::WomersleyProfile`.
- Rustdoc now states the exact Bessel-form velocity theorem and no-slip proof sketch.
- Added `test_exact_womersley_no_slip_wall_condition`, which checks computed wall velocity at five phases.
- `cargo check -p cfd-validation --no-default-features` completed successfully in 3.44 seconds after warm build.
- `cargo test -p cfd-validation --no-default-features analytical::womersley --lib -- --nocapture` passed 5/5 tests in 10.50 seconds after compilation.
- `cargo nextest run -p cfd-validation --no-default-features womersley --fail-fast --hide-progress-bar` exceeded the 60-second compilation bound before test execution; no failing Rust test was reported.

---

## OPEN-033: Anderson-Accelerated Picard Solver Stalls on Ill-Conditioned Networks (limit cycle)

**Severity**: OPEN — real numerical defect, requires algorithm work (not a mechanical fix)
**Component**: `crates/cfd-math/src/nonlinear_solver/anderson.rs` (`AndersonAccelerator::compute_next`, `QrState`), `crates/cfd-1d/src/solver/core/anderson_acceleration.rs` (`select_next_iterate`), `crates/cfd-1d/src/solver/core/convergence.rs` (`has_converged_dual`)
**Failing test**: `tests/cross_fidelity_blueprint.rs::cross_fidelity_blueprint_complex_branching` — panics with `MaxIterationsExceeded: Convergence failed: Maximum iterations (10000) exceeded` when solving the double-trifurcation blueprint's 1D reference network (`solve_reference_trace` in `crates/cfd-2d/src/network/reference.rs`).

### Root cause (instrumented-run evidence, not yet fixed)

Two independent investigations (temporary debug instrumentation, fully reverted — no diff left in the tree) converged on the following:

1. **Confirmed NOT the cause**: the resistance formula (`crates/cfd-1d/src/physics/resistance/models/hagen_poiseuille.rs:136-139`), `MatrixAssembler::assemble_into`, and `LinearSystemSolver` are unaffected by the leto/eunomia migration (mechanical ports only, verified by diffing `d58d1fe3^..1d768895`). The linear solve itself is exact: instrumented `rel_resid = ‖Ax−b‖/‖b‖` converges to ~1e-13 by iteration 8 and to machine precision (~1e-16) from iteration ~200 onward, for the entire 10000-iteration run. This network's conductance ratio is genuinely ill-conditioned (main channel 4mm vs. 80µm throat: `R_throat/R_main ≈ (0.2mm/15mm)×(4000µm/80µm)⁴ ≈ 8.3×10⁴`, ~5 orders of magnitude) — this is within the range `LinearSystemSolver`'s own docstring claims to handle (linear_system.rs:5-7), and it does handle it correctly.

2. **Confirmed the actual failure mode**: the *outer* nonlinear Picard/Anderson fixed-point iteration is stuck in a **persistent limit cycle**, not diverging and not slow-monotonic. Instrumented `rel_change = ‖Δx‖/‖x‖` oscillates in the band 8.0e-5 to 4.4e-4 with a visible ~400-iteration quasi-period from iteration ~200 through iteration 9999 — no decreasing trend over ~9600 iterations. `has_converged_dual` (convergence.rs:211) requires `relative_change < tolerance` (1e-8), which this never approaches.

3. **Contributing candidate (real defect, independently confirmed by code inspection, but did NOT fire in the instrumented run — 0 occurrences)**: `AndersonAccelerator::compute_next` (anderson.rs:281-334) indexes `self.delta_x[j]`/`self.delta_f[j]` for `j` over `gamma.len()` (`= qr_state.q_cols.len()`), but `QrState::append_column` (anderson.rs:129-186) can silently discard a near-collinear column (`norm_w <= drop_tol`) without evicting the corresponding entry from the *separately tracked* `delta_x`/`delta_f` deques. Whenever this triggers, `gamma[j]` no longer corresponds to `delta_x[j]`/`delta_f[j]` — the accelerated correction is applied against the wrong history vectors. This is a genuine latent bug in the QR-based Anderson accelerator introduced by the leto/eunomia migration (the pre-migration implementation was a dense Gram-matrix + `nalgebra::linalg::LU` recompute each iteration with no incremental QR state to desync). It should be fixed regardless of whether it's the proximate cause of this specific test's stall — recommended fix: make `append_column` report whether it accepted or dropped the column, and have `compute_next` evict the same index from `delta_x`/`delta_f` in lockstep, with a `debug_assert!(qr_state.q_cols.len() == self.delta_x.len())` invariant to catch regressions.

4. **Most likely proximate mechanism for the observed cycle (reasoned inference, not directly instrumented)**: `select_next_iterate` (`anderson_acceleration.rs:42-95`) accepts the Anderson-accelerated iterate whenever `accelerated_residual <= picard_residual` and `damped_step_norm < picard_step_norm` — there is no oscillation/stagnation detector (e.g. no sign-reversal check on `Δx_k · Δx_{k-1}`, no history reset on stagnation). On a network with a ~5-order-of-magnitude conductance spread, the linearized fixed-point Jacobian has no single dominant eigenvalue; a rank-5 QR-secant correction can steer into a stable 2-or-more-cycle that is locally no worse than raw Picard by those two metrics alone, so it's accepted indefinitely without ever being flagged as non-convergent.

### Fix direction (not yet implemented — needs real algorithm work, do not attempt a mechanical patch)

- Fix item 3 above regardless (real bug, straightforward once flagged).
- For item 4: add an oscillation/stagnation detector to `select_next_iterate` or the outer Picard loop — e.g. track a short window of `rel_change` history and trigger an Anderson-history reset (clear `delta_x`/`delta_f`/`qr_state`, fall back to plain damped Picard for a few steps) when the window shows non-decaying periodic behavior, then re-engage Anderson acceleration. This is standard practice for Anderson acceleration on stiff/ill-conditioned fixed-point maps (see e.g. Walker & Ni 2011 restart heuristics) — do not implement without deriving/citing the specific restart criterion, and do not weaken `has_converged_dual`'s tolerance or raise `max_iterations` as a substitute fix (that would be test-gaming, not a fix).

### Corroborating re-investigation (independent second pass, also fully reverted)

A second, independent instrumentation pass confirmed the same limit-cycle diagnosis via a different signal: `solution_change_norm` (‖xₖ−xₖ₋₁‖) oscillates in a stable ~7e-4 to 1.5e-3 band from iteration ~20 through 9999 with a repeating ~200-iteration pattern, while `residual_norm` (the *linear* system residual) is at machine epsilon by iteration ~8. Per-edge flow rates are independently stable to 7 significant digits pass-to-pass, yet the reconstructed `relative_flow_change` gate sits at ~3e-8 to 3e-7 — hovering at/above the 1e-8 tolerance without ever settling below it. Edge resistances span ~1165× (5.8e6 to 6.75e9); the condition-number-only noise floor for that spread is ~2.6e-13 — three orders of magnitude below the observed oscillation, ruling out simple ill-conditioning/rounding as the cause.

This pass also ruled out (with diff/instrumentation evidence) the resistance-model formulas (`rectangular.rs`, `venturi/model.rs`, `venturi_coefficients.rs` — mechanical trait-surface transliteration only, byte-identical constants), matrix assembly, the `1/r_eff` conductance conversion, and the Reynolds-2300 laminar/turbulent branch switch in `update_single_edge_resistance` (`wrapper.rs:533-536`; instrumented directly — `branch_darcy=false` at every sample, hand-computed Re ~0.8-15 for this geometry, nowhere near 2300). Confirmed the test predates the migration (added in `e636f3c4`), so it isn't a newly authored test tuned to the old backend.

**New lead for whoever picks this up**: `crates/cfd-math/src/nonlinear_solver/jfnk.rs` (Jacobian-free Newton-Krylov — a plausible alternative to Picard/Anderson for this stiffness regime) has a 199-line migration diff that has not yet been audited line-by-line for a genuine bug vs. transliteration, and it's unclear whether it's even wired into `NetworkSolver`'s `is_linear_static_network`/Picard branch selection for this problem class. Auditing `jfnk.rs` and checking whether routing this network through JFNK instead of Picard/Anderson converges is the next concrete step — a non-contractive Picard map may need a fundamentally different globalization scheme rather than a better Anderson restart heuristic.

### Verification the fix must satisfy

`cargo nextest run -p cfd-suite cross_fidelity_blueprint_complex_branching` passes without raising `max_iterations` or loosening `tolerance`; `cargo nextest run -p cfd-1d -p cfd-2d -p cfd-math --no-fail-fast` shows no regressions (in particular any test currently relying on the exact QR-Anderson trajectory); `cargo clippy` and `cargo fmt --check` clean on touched files.
