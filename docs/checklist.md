# Sprint 1.96.167 Checklist: cfd-math native sparse-LU result ownership
**Goal**: Consume the provider-owned native array view/result boundary from
the direct solver without redundant consumer-side allocations.

**Success Criteria**:
- [x] `direct_solver.rs` passes `rhs.view()` to `SparseLuSolver::solve_view`.
- [x] The primary path returns the provider-owned `Array1` directly.
- [x] Existing positive, singular, fallback, and generic scalar value tests
      pass without weakening assertions.
- [x] Direct-solver Rustdoc, changelog, and gap evidence describe the native
      ownership boundary and its evidence limits.
- [x] Focused format, check, warning-denied Clippy, Nextest, doctest, and
      Rustdoc gates pass.

**Residual**: provider branch/API integration and exact-head downstream
verification remain open until both repositories are delivered.

---

# Sprint 1.96.166 Checklist: cfd-math IncompleteCholesky Leto CSR
**Goal**: Move IncompleteCholesky construction and factors to Leto CSR.

**Success Criteria**:
- [x] `cholesky.rs` stores `leto_ops::CsrMatrix`.
- [x] Symmetry and factorization reads use Leto CSR row-value lookup.
- [x] IC(0) factors are constructed with `CsrMatrix::from_parts`.
- [x] Forward/backward substitution walks Leto CSR row views.
- [x] Source and integration preconditioner edge tests pass Leto CSR into
  `IncompleteCholesky::new`.
- [x] `cholesky.rs` and `tests/preconditioner_edge_cases.rs` have no direct
  nalgebra sparse/direct CSR/row-offset/get-entry residue for Cholesky.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run cfd-math all-target check.
- [x] Run lib/tests and all-target clippy.
- [x] Run cholesky-filter nextest (5/5 passed).
- [x] Run preconditioner-filter nextest (76/76 passed).
- [x] Run targeted Cholesky residue scans.

**Residual**: Schwarz, direct solver, remaining transitional solver fixtures,
and the shared `crate::sparse::SparseMatrix` solver matrix boundary still use
nalgebra-sparse.

---

# Sprint 1.96.165 Checklist: cfd-math ILU Leto CSR
**Goal**: Move ILU preconditioner family construction and factors to Leto CSR.

**Success Criteria**:
- [x] `preconditioners/ilu` stores `leto_ops::CsrMatrix`.
- [x] ILU(0) and ILU(k) factorization construct factors with Leto CSR
  `from_parts`.
- [x] ILU triangular solves walk Leto CSR row views.
- [x] Source tests, integration edge tests, `LinearSolverChain`, and Schwarz
  local solver setup pass Leto CSR into ILU or convert once at a remaining
  boundary.
- [x] `preconditioners/ilu` has no direct nalgebra sparse/direct
  CSR/row-offset/get-entry residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run cfd-math all-target check.
- [x] Run lib/tests and all-target clippy.
- [x] Run ilu-filter nextest (21/21 passed).
- [x] Run preconditioner-filter nextest (76/76 passed).
- [x] Run linear_solver::tests nextest (53/53 passed).
- [x] Run targeted ILU residue scans.

**Residual**: Schwarz, direct solver, remaining transitional solver fixtures,
and the shared `crate::sparse::SparseMatrix` solver matrix boundary still use
nalgebra-sparse.

---

# Sprint 1.96.164 Checklist: cfd-math SSOR Leto CSR
**Goal**: Move SSOR preconditioner construction to Leto CSR.

**Success Criteria**:
- [x] `ssor.rs` stores `leto_ops::CsrMatrix`.
- [x] SSOR forward/backward sweeps use Leto CSR row views.
- [x] Source preconditioner edge tests convert legacy solver CSR fixtures once
  before constructing SSOR.
- [x] `ssor.rs` has no direct nalgebra sparse/direct CSR/row-offset/get-entry
  residue.
- [x] SSOR vector-length regression asserts exact typed error messages.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run cfd-math all-target check.
- [x] Run lib/tests and all-target clippy.
- [x] Run ssor-filter nextest (5/5 passed).
- [x] Run preconditioner-filter nextest (76/76 passed).
- [x] Run targeted SSOR residue scans.

**Residual**: Schwarz, direct solver, integration-test fixtures,
and the shared `crate::sparse::SparseMatrix` solver matrix boundary still use
nalgebra-sparse.

---

# Sprint 1.96.163 Checklist: cfd-math Basic Preconditioner Leto CSR
**Goal**: Move Jacobi/SOR basic preconditioner construction to Leto CSR.

**Success Criteria**:
- [x] `basic.rs` constructs Jacobi from `leto_ops::CsrMatrix` diagonal access.
- [x] `basic.rs` constructs and applies SOR from `leto_ops::CsrMatrix` row
  access.
- [x] Source linear-solver tests and `core_solver_tests.rs` convert legacy
  solver CSR fixtures once before constructing Jacobi/SOR.
- [x] `basic.rs` has no direct nalgebra sparse/direct CSR/row-offset/get-entry
  residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run core solver test check.
- [x] Run cfd-math all-target check.
- [x] Run lib/tests, core solver, and all-target clippy.
- [x] Run linear_solver::tests nextest (53/53 passed).
- [x] Run core_solver_tests nextest (4/4 passed).
- [x] Run preconditioner-filter nextest (76/76 passed).
- [x] Run targeted basic-preconditioner residue scans.

**Residual**: Schwarz, direct solver, and the shared
`crate::sparse::SparseMatrix` solver matrix boundary still use nalgebra-sparse.

---

# Sprint 1.96.162 Checklist: cfd-math AMG/Coarsening Leto CSR
**Goal**: Move AMG/coarsening sparse storage and benchmark/test construction to
Leto CSR.

**Success Criteria**:
- [x] Multigrid AMG uses `leto_ops::CsrMatrix` for its local sparse boundary.
- [x] Coarsening/interpolation/smoothers/cycles use Leto CSR row access, SpMV,
  transpose, and sparse products.
- [x] `coarsening_bench.rs` and `algebraic_distance_bench.rs` construct Leto
  CSR directly.
- [x] AMG coarsening and integration tests construct AMG preconditioners from
  Leto CSR.
- [x] Migrated AMG/coarsening files have no direct nalgebra sparse/Coo/direct
  CSR/get_entry/old SpMV residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused AMG/coarsening test and bench checks.
- [x] Run focused clippy for cfd-math lib, touched tests, and touched benches.
- [x] Run AMG integration nextest (5/5 passed).
- [x] Run AMG-filter nextest (6/6 passed).
- [x] Run multigrid::coarsening nextest (10/10 passed).
- [x] Run targeted AMG/coarsening sparse residue scan.

**Residual**: the broader `crate::sparse::SparseMatrix`
solver/direct/preconditioner boundary still resolves to nalgebra-sparse.
`LinearSolverChain` converts that boundary once before AMG construction.

---

# Sprint 1.96.161 Checklist: cfd-math Leto CSR Benchmarks
**Goal**: Move SpMV/CG benchmark CSR storage off nalgebra sparse and onto Leto
CSR.

**Success Criteria**:
- [x] `spmv_bench.rs` constructs `leto_ops::CsrMatrix` directly.
- [x] `cg_bench.rs` constructs `leto_ops::CsrMatrix` directly.
- [x] `math_benchmarks.rs` CG benchmark constructs `leto_ops::CsrMatrix`
  directly.
- [x] SpMV benchmark dispatches through direct Leto CSR `LinearOperator::apply`.
- [x] Migrated benchmark files have no `nalgebra_sparse`, `CooMatrix`,
  `DVector`, `DMatrix`, `cfd_math::sparse::spmv`, or `try_from_csr_data`
  residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run focused checks for `spmv_bench`, `cg_bench`, and `math_benchmarks`.
- [x] Run focused clippy for all three migrated benches.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run sparse-filter nextest (19/19 passed).
- [x] Run targeted migrated-benchmark residue scan.

**Residual**: after Sprint 1.96.162 moved the AMG/coarsening boundary to Leto
CSR, the remaining sparse provider gap is the broader
solver/direct/preconditioner nalgebra-sparse matrix surface.

---

# Sprint 1.96.160 Checklist: cfd-math Leto CSR LinearOperator
**Goal**: Add a direct Leto CSR sparse-operator path for cfd-math solvers.

**Success Criteria**:
- [x] `leto_ops::CsrMatrix<T>` implements cfd-math `LinearOperator<T>`.
- [x] Legacy nalgebra CSR SpMV delegates through the same Leto CSR helper.
- [x] Simple GMRES integration tests build Leto CSR matrices directly.
- [x] Simple GMRES residual checks use the direct Leto CSR operator path.
- [x] `simple_gmres_tests.rs` has no `nalgebra_sparse`, `CooMatrix`,
  `nalgebra::`, `DVector`, `DMatrix`, or `cfd_math::sparse` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run simple_gmres test check.
- [x] Run focused simple_gmres nextest (3/3 passed).
- [x] Run focused simple_gmres clippy.
- [x] Run cfd-math lib check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run sparse-filter nextest (19/19 passed).
- [x] Run gmres-filter nextest (21/21 passed).
- [x] Run targeted simple GMRES residue scan.

**Residual**: cfd-math still exposes the broader nalgebra-sparse matrix
boundary in sparse builders, preconditioners, AMG, direct solver, tests, and
benches. This slice opens and verifies the direct Leto CSR solver path.

---

# Sprint 1.96.159 Checklist: cfd-math Storage-Slice Closure
**Goal**: Close the remaining cfd-math Leto storage-slice residue.

**Success Criteria**:
- [x] Nonlinear dense pivoting uses direct `Array1`/`Array2` indexing.
- [x] Nonlinear linalg exposes no `StorageMut` mutable slice helpers.
- [x] AMG interpolation quality validation indexes interpolated vectors
  directly.
- [x] Multigrid smoother tests assert values through direct indexing.
- [x] cfd-math `src`/`tests` has no `leto::Storage`, `StorageMut`,
  `.storage().as_slice()`, `as_slice_mut()`, `vector_slice_mut`, or
  `matrix_slice_mut` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused cfd-math nonlinear_solver/multigrid nextest (46/46 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run cfd-math source/test storage-slice residue scan.

**Residual**: cfd-math Leto storage-slice cleanup is closed. Remaining
cfd-math Atlas-provider work is nalgebra/nalgebra-sparse replacement and other
provider boundaries.

---

# Sprint 1.96.158 Checklist: cfd-math Sparse/Basic Leto Array1
**Goal**: Remove direct Leto storage-slice dependencies from sparse operations
and the basic Jacobi preconditioner.

**Success Criteria**:
- [x] SpMV input access uses direct `Array1` indexing.
- [x] SpMV output writeback uses direct `Array1` indexing.
- [x] Row scaling uses direct `Array1` indexing.
- [x] Column scaling uses direct `Array1` indexing.
- [x] Jacobi preconditioner diagonal iteration uses direct `Array1` indexing.
- [x] Migrated sparse/basic files have no `leto::Storage`,
  `.storage().as_slice()`, `as_slice_mut()`, or obsolete SpMV output
  contiguity diagnostic residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused cfd-math sparse/preconditioner nextest (95/95 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run targeted sparse/basic residue scans.

**Residual**: remaining cfd-math storage-slice owners are nonlinear mutable
dense-workspace helpers and multigrid interpolation/smoother internals.

---

# Sprint 1.96.157 Checklist: cfd-math GPU Operator Leto Array1
**Goal**: Remove the GPU linear operator's direct dependency on Leto storage
slices while preserving the Hephaestus-backed GPU execution path.

**Success Criteria**:
- [x] GPU operator input access uses direct `Array1` indexing.
- [x] GPU operator output writeback uses direct `Array1` indexing.
- [x] GPU operator rejects mismatched input/output vector lengths before GPU
  upload.
- [x] `operators/gpu.rs` has no `leto::Storage`,
  `.storage().as_slice()`, `as_slice_mut()`, or obsolete output-contiguity
  diagnostic residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math GPU-feature check.
- [x] Run focused cfd-math GPU-feature linear_solver::operators nextest (5/5
  passed).
- [x] Run cfd-math GPU-feature lib clippy.
- [x] Run cfd-math GPU-feature all-target check.
- [x] Run cfd-math GPU-feature all-target clippy.
- [x] Run targeted GPU operator residue scan.

**Residual**: remaining cfd-math storage-slice owners are sparse operations and
multigrid internals. Broader GPU provider work remains in raw WGPU
context/kernel ownership outside this operator.

---

# Sprint 1.96.156 Checklist: cfd-math Finite-Difference Operators Leto Array1
**Goal**: Remove Leto storage-slice dependencies from CPU finite-difference
linear operators.

**Success Criteria**:
- [x] 2D Laplacian input/output access uses direct `Array1` indexing.
- [x] 3D Poisson input/output access uses direct `Array1` indexing.
- [x] 1D momentum input/output access uses direct `Array1` indexing.
- [x] 2D momentum input/output access uses direct `Array1` indexing.
- [x] 2D energy input/output access uses direct `Array1` indexing.
- [x] Migrated operator files have no `leto::Storage`,
  `.storage().as_slice()`, or `as_slice_mut()` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused cfd-math linear_solver::operators nextest (5/5 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run targeted finite-difference operator residue scan.

**Residual**: remaining storage-slice owners are sparse operations, GPU
operator, and multigrid internals.

---

# Sprint 1.96.155 Checklist: cfd-math Nonlinear Linalg Leto Array1
**Goal**: Remove the nonlinear solver helper module's immutable dependency on
Leto storage slices.

**Success Criteria**:
- [x] `dot` uses direct `Array1` indexing.
- [x] `add` uses direct `Array1` indexing.
- [x] `sub` uses direct `Array1` indexing.
- [x] `add_scaled` uses direct `Array1` indexing.
- [x] `add_scaled_in_place` reads RHS through direct `Array1` indexing.
- [x] `scale` uses direct `Array1` indexing.
- [x] Anderson acceleration indexes the `gamma` vector directly.
- [x] `linalg.rs` has no immutable `.storage().as_slice()` or `leto::Storage`
  residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused cfd-math nonlinear_solver nextest (9/9 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run targeted nonlinear linalg/Anderson residue scan.

**Residual**: nonlinear linalg mutable dense-workspace helpers still use
`StorageMut`; remaining immutable storage-slice owners are sparse operations,
linear operators, GPU operator, and multigrid internals.

---

# Sprint 1.96.154 Checklist: cfd-math Production SIMD Vector Leto Array1
**Goal**: Remove the production SIMD vector extension's direct dependency on
Leto storage slices while preserving Moirai execution dispatch.

**Success Criteria**:
- [x] `simd_mul` uses direct `Array1` indexing.
- [x] `simd_dot` uses direct `Array1` indexing.
- [x] `simd_norm` uses direct `Array1` indexing.
- [x] `par_map` uses direct `Array1` indexing.
- [x] `src/simd/vector.rs` has no `leto::Storage` import.
- [x] `src/simd/vector.rs` has no `.storage().as_slice()` calls.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run focused cfd-math simd::vector nextest (1/1 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run broader cfd-math SIMD-filter nextest (26/26 passed).
- [x] Run targeted source-residue scan.

**Residual**: source-level storage-slice owners remain in nonlinear linalg,
sparse operations, linear operators, GPU operator, and multigrid internals.

---

# Sprint 1.96.153 Checklist: cfd-math SIMD Integration Test Leto Array1
**Goal**: Remove the remaining cfd-math integration-test Leto storage-slice
bridge from `simd_tests.rs`.

**Success Criteria**:
- [x] `test_cfd_residual_calculation` reads the Leto `spmv` result through
  direct `Array1` indexing.
- [x] `simd_tests.rs` has no `leto::Storage` import.
- [x] `simd_tests.rs` has no `.storage().as_slice()` conversion.
- [x] `crates/cfd-math/tests` has no `DVector`, nalgebra vector import, local
  preconditioner bridge, `Storage`, or storage-slice conversion residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run simd_tests check.
- [x] Run focused cfd-math simd_tests nextest (12/12 passed).
- [x] Run simd_tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run broader cfd-math SIMD-filter nextest (26/26 passed).
- [x] Run targeted cfd-math integration-test provider-residue scans.

**Residual**: remaining provider work is outside the cfd-math integration-test
vector bridge layer: `nalgebra_sparse::CsrMatrix`, dense nalgebra test oracles,
and source-level Leto storage-slice internals.

---

# Sprint 1.96.152 Checklist: cfd-math AMG Integration Test Leto Array1
**Goal**: Move cfd-math AMG integration vector paths off nalgebra `DVector`
bridge helpers and onto direct Leto array inputs.

**Success Criteria**:
- [x] Exact solution helper returns `leto::Array1`.
- [x] RHS assembly returns `leto::Array1`.
- [x] BiCGSTAB AMG integration test solves and checks error on Leto arrays.
- [x] GMRES AMG integration test solves and checks error on Leto arrays.
- [x] AMG cycle comparison uses direct Leto preconditioner application.
- [x] Two-grid convergence factor applies AMG to Leto work buffers.
- [x] `amg_integration_test.rs` has no `DVector`, nalgebra vector import,
  `Storage`, storage-slice conversion, or local preconditioner bridge residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run amg_integration_test check.
- [x] Run focused cfd-math amg_integration_test nextest (5/5 passed).
- [x] Run amg_integration_test clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run broader cfd-math AMG-filter nextest (6/6 passed).
- [x] Run targeted amg_integration_test provider-residue scan.

**Residual**: `amg_integration_test.rs` still uses
`nalgebra_sparse::CsrMatrix` for sparse storage and nalgebra
`DMatrix`/`SymmetricEigen` for the dense energy-norm oracle. Remaining
cfd-math integration storage-slice residue is in `tests/simd_tests.rs`.

---

# Sprint 1.96.151 Checklist: cfd-math Preconditioner Edge-Case Tests Leto Array1
**Goal**: Move cfd-math preconditioner edge-case integration tests off
nalgebra `DVector` bridge helpers and onto direct Leto array inputs.

**Success Criteria**:
- [x] ILU(0) ill-conditioned test uses `leto::Array1`.
- [x] ILU(k) conditioning test uses `leto::Array1`.
- [x] Repeated ILU application test uses `leto::Array1`.
- [x] Extreme-value ILU test uses `leto::Array1`.
- [x] ILU sparsity-preservation test uses `leto::Array1`.
- [x] `preconditioner_edge_cases.rs` has no `DVector`, nalgebra vector import,
  `Storage`, storage-slice conversion, or local preconditioner bridge residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run preconditioner_edge_cases check.
- [x] Run focused cfd-math preconditioner_edge_cases nextest (6/6 passed).
- [x] Run preconditioner_edge_cases clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run broader cfd-math preconditioner nextest (76/76 passed).
- [x] Run targeted preconditioner_edge_cases provider-residue scan.

**Residual**: `preconditioner_edge_cases.rs` still uses the shared
`nalgebra_sparse::CsrMatrix` matrix boundary. Remaining cfd-math integration
storage-slice residue is in `tests/simd_tests.rs`.

---

# Sprint 1.96.150 Checklist: cfd-math Linear-Solver Test Module Leto Array1
**Goal**: Move the cfd-math linear-solver source test module tree off
nalgebra `DVector` bridge helpers and onto direct Leto array inputs.

**Success Criteria**:
- [x] `src/linear_solver/tests/mod.rs` uses `leto::Array1`.
- [x] `edge_case_tests.rs` uses `leto::Array1`.
- [x] `adversarial_solver_tests.rs` uses `leto::Array1`.
- [x] `extended_edge_case_tests.rs` uses `leto::Array1`.
- [x] Residual checks use the Leto SpMV/helper path.
- [x] The searched solver/sparse cone has no `T::zero()`, `T::one()`,
  `default_epsilon()`, or `to_subset()` residue.
- [x] `src/linear_solver/tests` has no `DVector`, nalgebra vector import,
  `Storage`, storage-slice conversion, local solve/preconditioner bridge,
  matrix-vector `&a * &x`, or vector `.norm()` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math lib check.
- [x] Run cfd-math all-target check.
- [x] Run focused cfd-math linear_solver::tests nextest (53/53 passed).
- [x] Run broader cfd-math linear_solver nextest (176/176 passed).
- [x] Run cfd-math lib/tests clippy.
- [x] Run cfd-math all-target clippy.
- [x] Run targeted source-test provider-residue scan.
- [x] Run targeted old scalar helper residue scan.

**Residual**: integration-test vector bridge residue was closed by Sprint
1.96.153; production sparse/scalar boundaries still include
`nalgebra_sparse::CsrMatrix` and transitional nalgebra scalar/test utilities.

---

# Sprint 1.96.149 Checklist: cfd-math Core Solver Tests Leto Array1
**Goal**: Move cfd-math core solver validation tests off nalgebra `DVector`
bridge helpers and onto direct Leto array inputs.

**Success Criteria**:
- [x] BiCGSTAB validation test uses `leto::Array1`.
- [x] GMRES validation test uses `leto::Array1`.
- [x] Preconditioner integration test uses `leto::Array1`.
- [x] Condition-number robustness test uses `leto::Array1`.
- [x] Residual checks use the Leto SpMV path with explicit thresholds.
- [x] `core_solver_tests.rs` has no `DVector`, nalgebra vector, `Storage`,
  storage-slice conversion, local solve/preconditioner bridge, matrix-vector
  `&a * &x`, or vector `.norm()` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run core_solver_tests check.
- [x] Run focused cfd-math core_solver_tests nextest (4/4 passed).
- [x] Run core_solver_tests clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run targeted core solver test provider-residue scan.

**Residual**: `core_solver_tests.rs` still uses the shared
`nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix boundary, and broader cfd-math
integration/adversarial/preconditioner tests still contain nalgebra `DVector`
bridge diagnostics.

---

# Sprint 1.96.148 Checklist: cfd-math Simple GMRES Tests Leto Array1
**Goal**: Move cfd-math simple GMRES integration tests off nalgebra `DVector`
bridge macros and onto direct Leto array inputs.

**Success Criteria**:
- [x] Basic GMRES integration test uses `leto::Array1`.
- [x] Restarted GMRES integration test uses `leto::Array1`.
- [x] Preconditioned GMRES integration test uses `leto::Array1`.
- [x] Residual checks use the Leto SpMV path with explicit thresholds.
- [x] `simple_gmres_tests.rs` has no `DVector`, nalgebra vector,
  `Storage`, storage-slice conversion, local solve-bridge macro,
  matrix-vector `&a * &x`, or vector `.norm()` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run simple_gmres test check.
- [x] Run focused cfd-math simple_gmres nextest (3/3 passed).
- [x] Run simple_gmres clippy.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run touched-file `git diff --check`.
- [x] Run targeted simple GMRES test provider-residue scan.

**Residual**: `simple_gmres_tests.rs` still uses the shared
`nalgebra_sparse::CsrMatrix`/`CooMatrix` matrix boundary, and broader cfd-math
integration/adversarial/core/preconditioner tests still contain nalgebra
`DVector` bridge diagnostics.

---

# Sprint 1.96.147 Checklist: cfd-math Matrix-Free Tests Leto Array1
**Goal**: Move cfd-math matrix-free solver tests off nalgebra `DVector` bridge
macros and onto direct Leto array inputs.

**Success Criteria**:
- [x] Matrix-free CG identity test uses `leto::Array1`.
- [x] Matrix-free GMRES identity test uses `leto::Array1`.
- [x] Scaled-operator integration test uses `leto::Array1`.
- [x] Operator-size mismatch test asserts the exact typed
  `InvalidConfiguration` message.
- [x] `matrix_free/tests.rs` has no `DVector`, nalgebra, `Storage`,
  storage-slice conversion, or local solve-bridge macro residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math matrix-free nextest (4/4 passed).
- [x] Run touched-file `git diff --check`.
- [x] Run targeted matrix-free test provider-residue scan.

**Residual**: broader cfd-math linear-solver integration/adversarial/core/
preconditioner test diagnostics still contain nalgebra `DVector` bridges, and
production sparse/scalar provider boundaries still include
`nalgebra_sparse::CsrMatrix` and transitional `nalgebra::RealField`.

---

# Sprint 1.96.146 Checklist: cfd-math BiCGSTAB Leto Array1
**Goal**: Move cfd-math BiCGSTAB solver workspaces, direct solve methods, and
the final solver-chain fallback off nalgebra `DVector` conversion and onto Leto
arrays.

**Success Criteria**:
- [x] BiCGSTAB residual/search/operator workspaces use `leto::Array1`.
- [x] `solve_preconditioned` accepts Leto RHS/solution arrays.
- [x] `solve_unpreconditioned` accepts Leto RHS/solution arrays.
- [x] `IterativeLinearSolver::solve` dispatches without legacy vector bridge
  helpers.
- [x] `LinearSolverChain` final BiCGSTAB fallback stays on Leto arrays.
- [x] CG and BiCGSTAB share one Leto vector helper module.
- [x] Obsolete legacy linear-solver vector bridge helpers are removed.
- [x] Migrated BiCGSTAB/chain/traits source has no `DVector`, legacy bridge
  helper, nalgebra vector math, or nalgebra `copy_from` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math BiCGSTAB nextest (24/24 passed).
- [x] Run broader cfd-math linear-solver nextest (176/176 passed).
- [x] Run cfd-math AMG integration nextest (5/5 passed).
- [x] Run touched-file `git diff --check`.
- [x] Run targeted BiCGSTAB/chain/traits provider-residue scans.

**Residual**: cfd-math still carries the shared `nalgebra_sparse::CsrMatrix`
storage/provider boundary, transitional `nalgebra::RealField` scalar bounds,
and nalgebra `DVector` in remaining matrix-free/preconditioner/integration test
diagnostics.

---

# Sprint 1.96.145 Checklist: cfd-math Conjugate Gradient Leto Array1
**Goal**: Move cfd-math Conjugate Gradient solver workspaces and direct solve
methods off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] CG reusable workspaces use `leto::Array1`.
- [x] `solve_preconditioned` accepts Leto RHS/solution arrays.
- [x] `solve_unpreconditioned` accepts Leto RHS/solution arrays.
- [x] `IterativeLinearSolver::solve` dispatches without legacy vector bridge
  helpers.
- [x] CG benchmark call sites construct Leto vectors at the measured API.
- [x] Co-located tests assert solved values and exact dimension/max-iteration
  errors.
- [x] Migrated CG code and CG benchmark call sites have no `DVector`,
  legacy bridge helper calls, `num_traits`, storage-slice conversion, nalgebra
  vector math calls, or nalgebra `copy_from` residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math conjugate nextest (13/13 passed).
- [x] Run broader cfd-math linear-solver nextest (176/176 passed).
- [x] Run targeted CG provider-residue scan.

**Residual after the BiCGSTAB follow-up**: the shared linear-solver trait
family still carries the transitional `nalgebra::RealField` scalar bound and
sparse storage remains on `nalgebra_sparse::CsrMatrix`.

---

# Sprint 1.96.144 Checklist: cfd-math Schwarz Preconditioner Leto Array1
**Goal**: Move cfd-math Schwarz local apply paths and local RHS workspaces off
nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] Additive Schwarz accepts and returns Leto vectors.
- [x] Multiplicative Schwarz accepts and returns Leto vectors.
- [x] Local RHS extraction uses Leto arrays directly.
- [x] `Preconditioner::apply_to` no longer constructs nalgebra `DVector`
  workspaces.
- [x] Residual/output length mismatches return typed configuration errors.
- [x] `schwarz.rs` has no `DVector`, `Storage`, `num_traits`,
  `FromPrimitive`, local-RHS conversion bridge, `DMatrix`, or `ndarray`
  residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math Schwarz nextest (3/3 passed).
- [x] Run broader cfd-math preconditioner nextest (76/76 passed).
- [x] Run targeted Schwarz provider-residue scan.

**Residual**: Schwarz still stores and constructs sparse matrices through the
shared `nalgebra_sparse::CsrMatrix` boundary and remains constrained by the
global `Preconditioner<T>` nalgebra scalar bound.

---

# Sprint 1.96.143 Checklist: cfd-math ILU Triangular Solve Leto Array1
**Goal**: Move cfd-math ILU residual, intermediate, and solution workspaces
off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] ILU forward substitution accepts Leto residual/intermediate buffers.
- [x] ILU backward substitution accepts Leto intermediate/output buffers.
- [x] `IncompleteLU::apply_to` no longer constructs nalgebra `DVector`
  workspaces.
- [x] Residual/output length mismatches return typed configuration errors.
- [x] U-solve diagonal identity is owned by Eunomia `NumericElement`.
- [x] `ilu/types.rs` and `ilu/triangular.rs` have no `DVector`, `DMatrix`,
  `Storage`, old scalar identities, fallback wording, `component_mul`,
  `rows_mut`, row-view residue, `num_traits`, or `num_complex`.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math ILU nextest (21/21 passed).
- [x] Run broader cfd-math preconditioner nextest (74/74 passed).
- [x] Run targeted ILU provider-residue scan.

**Residual**: IncompleteLU still stores the shared
`nalgebra_sparse::CsrMatrix` LU factor boundary and remains constrained by the
global `Preconditioner<T>` nalgebra scalar bound.

---

# Sprint 1.96.142 Checklist: cfd-math Deflation Preconditioner Leto Array1
**Goal**: Move cfd-math deflation eigenvector state and projection work off
nalgebra `DVector` and onto Leto arrays.

**Success Criteria**:
- [x] Deflation eigenvector storage uses `leto::Array1`.
- [x] `add_eigenpair` accepts Leto eigenvectors.
- [x] `DeflationPreconditioner::apply_to` no longer constructs nalgebra
  `DVector` workspaces.
- [x] Output/eigenvector length mismatches return typed configuration errors.
- [x] Zero eigenvalues are rejected before projection division.
- [x] `deflation.rs` has no `DVector`, `DMatrix`, `Storage`, old scalar
  identities, fallback wording, `component_mul`, `rows_mut`, row-view residue,
  `num_traits`, or `num_complex`.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math deflation nextest (3/3 passed).
- [x] Run broader cfd-math preconditioner nextest (73/73 passed).
- [x] Run targeted deflation provider-residue scan.

**Residual**: Deflation still wraps the base preconditioner behind the existing
`Box<dyn Preconditioner<T>>`, and the global `Preconditioner<T>` trait still
uses the transitional nalgebra `RealField` scalar bound.

---

# Sprint 1.96.141 Checklist: cfd-math Basic Preconditioners Leto Array1
**Goal**: Move cfd-math Identity/Jacobi/SOR basic preconditioner vector state
and apply paths off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] Identity applies directly from Leto residual to Leto output.
- [x] Jacobi inverse diagonal storage uses `leto::Array1`.
- [x] Jacobi applies directly from Leto residual to Leto output.
- [x] SOR applies directly from Leto residual to Leto output.
- [x] Identity/Jacobi/SOR length mismatches return typed configuration errors
  with value-semantic regression coverage.
- [x] Basic preconditioner scalar identities and diagonal tolerance route
  through Eunomia.
- [x] `basic.rs` has no `DVector`, `DMatrix`, old scalar identities, fallback
  wording, `component_mul`, `rows_mut`, or row-view residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math basic mismatch nextest (1/1 passed).
- [x] Run broader cfd-math preconditioner nextest (70/70 passed).
- [x] Run targeted basic-preconditioner provider-residue scan.

**Residual**: Jacobi and SOR still store the shared
`nalgebra_sparse::CsrMatrix` boundary and remain constrained by the global
`Preconditioner<T>` nalgebra scalar bound.

---

# Sprint 1.96.140 Checklist: cfd-math IncompleteCholesky Leto Array1
**Goal**: Move cfd-math IncompleteCholesky residual, intermediate, and
solution workspaces off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] Cholesky forward substitution accepts Leto residual/intermediate buffers.
- [x] Cholesky backward substitution accepts Leto intermediate/output buffers.
- [x] `IncompleteCholesky::apply_to` no longer constructs nalgebra `DVector`
  workspaces.
- [x] Residual/output length mismatches return typed configuration errors.
- [x] IC(0) square-root dispatch is owned by Eunomia `NumericElement`.
- [x] `cholesky.rs` has no `DVector`, `DMatrix`, `Storage`, old scalar
  identities, silent fallback wording, ambiguous `.sqrt()`, or residual
  `DVector` workspace construction.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math Cholesky nextest (5/5 passed).
- [x] Run targeted Cholesky provider-residue scan.

**Follow-up**: Sprint 1.96.166 moved IncompleteCholesky factor storage to Leto
CSR. Residual sparse-provider work is now Schwarz, direct solver, and the
shared solver matrix boundary.

---

# Sprint 1.96.139 Checklist: cfd-math SSOR Preconditioner Leto Array1
**Goal**: Move cfd-math SSOR sweep residual and solution workspaces off
nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] SSOR forward sweep accepts Leto residual/output buffers.
- [x] SSOR backward sweep accepts Leto residual/output buffers.
- [x] `SSOR::apply_to` no longer constructs nalgebra `DVector` workspaces.
- [x] Residual/output length mismatches return typed configuration errors.
- [x] `ssor.rs` has no `DVector`, `DMatrix`, `Storage`, nalgebra row-view
  operations, `component_mul`, old scalar identities, silent fallback wording,
  or clone fallback residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math SSOR nextest (5/5 passed).
- [x] Run targeted SSOR provider-residue scan.

**Residual**: SSOR still stores the shared `nalgebra_sparse::CsrMatrix`
matrix boundary and remains constrained by the global `Preconditioner<T>`
nalgebra scalar bound.

---

# Sprint 1.96.138 Checklist: cfd-math Block/SIMPLE Preconditioner Leto Array1
**Goal**: Move cfd-math block/SIMPLE preconditioner diagonal and Schur vector
state off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] Diagonal inverse storage uses `leto::Array1`.
- [x] SIMPLE Schur diagonal inverse storage uses `leto::Array1`.
- [x] Direct block/SIMPLE `apply` methods accept and return Leto arrays.
- [x] `Preconditioner::apply_to` implementations no longer construct
  nalgebra `DVector` bridges.
- [x] Mismatched vector lengths return typed configuration errors.
- [x] `block_preconditioner.rs` has no `DVector`, `DMatrix`, nalgebra import,
  nalgebra row-view operations, `component_mul`, old scalar identities, silent
  fallback wording, or clone fallback residue.

### Closure
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math block-preconditioner nextest (4/4 passed).
- [x] Run targeted block-preconditioner provider-residue scan.

**Residual**: the global `Preconditioner` trait still carries the transitional
`nalgebra::RealField` scalar bound, and sparse storage still aliases
`nalgebra_sparse::CsrMatrix`.

---

# Sprint 1.96.137 Checklist: cfd-math GMRES Leto Workspace
**Goal**: Move cfd-math GMRES work/residual vectors and solver-chain GMRES
tiers off nalgebra `DVector` conversion and onto Leto arrays.

**Success Criteria**:
- [x] GMRES Arnoldi basis-column extraction uses Leto workspace arrays.
- [x] GMRES work, preconditioned work, and residual checks use `leto::Array1`.
- [x] GMRES direct solve methods accept Leto RHS/solution arrays.
- [x] `LinearSolverChain` keeps GMRES+AMG, GMRES+block, unpreconditioned
  GMRES, and GMRES+ILU tiers on Leto arrays.
- [x] `linear_solver/gmres` has no `DVector` or legacy operator/preconditioner
  bridge residue.

### Closure
- [x] Run cfd-math lib/tests/all-targets checks.
- [x] Run cfd-math fmt check.
- [x] Run cfd-math all-target clippy.
- [x] Run focused cfd-math GMRES nextest (21/21 passed).
- [x] Run cfd-2d no-default lib check.
- [x] Run focused cfd-2d momentum nextest (53/53 passed).
- [x] Run targeted GMRES `DVector`/legacy bridge residue scan.

**Residual**: cfd-math still has nalgebra scalar bounds, CG/BiCGSTAB
workspaces, `LinearSolverChain`'s final BiCGSTAB bridge, and
`nalgebra_sparse::CsrMatrix` storage pending later Leto/Eunomia slices.

---

# Sprint 1.96.136 Checklist: cfd-2d Momentum Leto Array1
**Goal**: Move cfd-2d momentum RHS/solution vectors and scalar seam off direct
nalgebra ownership.

**Success Criteria**:
- [x] Momentum solver stores RHS and solution vectors as Leto arrays.
- [x] Momentum boundary helpers mutate Leto RHS arrays directly.
- [x] Momentum solve paths call Leto-native iterative and direct solver APIs.
- [x] `cfd-2d` has no direct source/manifest `nalgebra`,
  `nalgebra-sparse`, `DVector`, `DMatrix`, or `linear_solver_bridge` residue.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check.
- [x] Run cfd-2d no-default all-target clippy.
- [x] Run focused cfd-2d momentum nextest (53/53 passed).
- [x] Run direct cfd-2d source/manifest nalgebra residue scan.
- [x] Run cfd-2d nalgebra/nalgebra-sparse cargo-tree audit.

**Residual**: cfd-2d still resolves nalgebra/nalgebra-sparse transitively
through upstream crates including cfd-1d, cfd-core, cfd-math, cfd-schematics,
and Gaia.

---

# Sprint 1.96.135 Checklist: cfd-2d Pressure-Velocity Leto Array1
**Goal**: Move cfd-2d pressure-velocity pressure-correction vectors from
nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] Pressure-velocity caches store RHS and correction solution as Leto arrays.
- [x] Pressure-velocity dispatch calls Leto-native iterative and direct solver
  APIs directly.
- [x] Targeted scans find no `DVector`, nalgebra, or
  `linear_solver_bridge` residue under `crates/cfd-2d/src/pressure_velocity`.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check.
- [x] Run cfd-2d no-default all-target clippy.
- [x] Run focused cfd-2d pressure-velocity nextest (16/16 passed).
- [x] Run targeted pressure-velocity DVector/nalgebra/bridge residue scan.

**Residual**: direct cfd-2d source/manifest nalgebra ownership is now removed;
nalgebra remains transitive through upstream owners.

---

# Sprint 1.96.134 Checklist: cfd-2d SIMPLE Leto Array1
**Goal**: Move cfd-2d SIMPLE pressure-correction vectors from nalgebra
`DVector` to Leto `Array1`.

**Success Criteria**:
- [x] SIMPLE stores RHS and `p_prime` as Leto arrays.
- [x] SIMPLE pressure solve calls `IterativeLinearSolver::solve` directly.
- [x] Targeted scans find no `DVector`, nalgebra, or
  `linear_solver_bridge` residue under `crates/cfd-2d/src/solvers/simple`.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check.
- [x] Run cfd-2d no-default all-target clippy.
- [x] Run focused cfd-2d SIMPLE nextest (19/19 passed).
- [x] Run targeted SIMPLE DVector/nalgebra/bridge residue scan.

**Residual**: direct cfd-2d source/manifest nalgebra ownership is now removed;
nalgebra remains transitive through upstream owners.

---

# Sprint 1.96.133 Checklist: cfd-2d FDM Leto Array1
**Goal**: Move cfd-2d FDM RHS and Gauss-Seidel solution vectors from nalgebra
`DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `solve_gauss_seidel` accepts a Leto RHS and returns a Leto solution.
- [x] Poisson FDM stencil assembly mutates a Leto RHS array.
- [x] Advection-diffusion FDM stencil assembly mutates a Leto RHS array.
- [x] Targeted scans find no `DVector` or nalgebra residue under
  `crates/cfd-2d/src/solvers/fdm`.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check.
- [x] Run cfd-2d no-default all-target clippy.
- [x] Run focused cfd-2d FDM nextest (2/2 passed).
- [x] Run targeted FDM DVector/nalgebra residue scan.

**Residual**: direct cfd-2d source/manifest nalgebra ownership is now removed;
nalgebra remains transitive through upstream owners.

---

# Sprint 1.96.132 Checklist: cfd-2d Time Integration Leto Array1
**Goal**: Move cfd-2d time-integration state vectors from nalgebra `DVector`
to Leto `Array1`.

**Success Criteria**:
- [x] Explicit time schemes accept and return Leto state vectors.
- [x] Implicit and multistep schemes use Leto vector arithmetic and
  `leto_ops::norm_l2` convergence checks.
- [x] Adaptive controller and adaptive integrator expose Leto state vectors.
- [x] Time tests construct and assert Leto `Array1` states.
- [x] Targeted scans find no `DVector` or nalgebra residue under
  `crates/cfd-2d/src/schemes/time`.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check.
- [x] Run cfd-2d no-default all-target clippy.
- [x] Run focused cfd-2d time nextest (29/29 passed).
- [x] Run targeted time-cone DVector/nalgebra residue scans.

**Residual**: direct cfd-2d source/manifest nalgebra ownership is now removed;
nalgebra remains transitive through upstream owners.

---

# Sprint 1.96.131 Checklist: cfd-2d DMatrix Leto Array2
**Goal**: Move the compact cfd-2d `DMatrix` residue set to Leto `Array2`.

**Success Criteria**:
- [x] Immersed-boundary force/velocity matrix APIs use `leto::Array2`.
- [x] `schemes::Grid2D` stores field values in `leto::Array2`.
- [x] Scheme callers/tests use Leto `shape()` and `[[i, j]]` indexing.
- [x] `blood_venturi` example passes Leto arrays through IBM coupling.
- [x] Targeted scans find no cfd-2d `DMatrix` residue or nalgebra-style
  `Grid2D.data` tuple access.

### Closure
- [x] Run cfd-2d fmt check.
- [x] Run cfd-2d no-default lib check/clippy.
- [x] Run cfd-2d no-default `blood_venturi` example check.
- [x] Run focused cfd-2d nextest (60/60 passed).
- [x] Run targeted DMatrix and grid-access residue scans.

**Residual**: cfd-2d still retains nalgebra in non-DMatrix `DVector`/sparse
linear-system boundaries and other provider seams.

---

# Sprint 1.96.130 Checklist: cfd-1d Vascular Eunomia Complex
**Goal**: Move cfd-1d Bessel/Womersley complex arithmetic from nalgebra
`Complex`/`ComplexField` to Eunomia's complex scalar provider.

**Success Criteria**:
- [x] Bessel recurrence accepts and returns `eunomia::Complex`.
- [x] Womersley profile constructs Eunomia complex values for analytical
  velocity, wall-shear, and flow-rate evaluation.
- [x] Convergence uses Eunomia complex `norm()` instead of
  `ComplexField::modulus()`.
- [x] Targeted vascular scans find no nalgebra complex residue.

### Closure
- [x] Run cfd-1d lib check/clippy.
- [x] Run focused Bessel/Womersley nextest (26/26 passed).
- [x] Run cfd-1d fmt check.
- [x] Run targeted nalgebra-complex residue scan.

**Residual**: cfd-1d still depends on nalgebra for network sparse/dense
linear-system storage and the transitional scalar seam.

---

# Sprint 1.96.129 Checklist: LinearOperator Trait Leto Vectors
**Goal**: Move the public `LinearOperator::apply` vector boundary from
nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `LinearOperator::apply` accepts Leto input and output arrays.
- [x] `LinearOperator::apply_transpose` follows the same Leto boundary.
- [x] cfd-math operator implementations compile through the Leto public
  boundary.
- [x] Solver internals bridge only at current nalgebra workspaces.
- [x] Targeted scans find no old public operator DVector signature.

### Closure
- [x] Run cfd-math all-target check/clippy.
- [x] Run focused cfd-math solver/operator nextest (80/80 passed).
- [x] Run cfd-math fmt check.
- [x] Run targeted DVector operator-signature residue scan.

**Residual after the CG/BiCGSTAB follow-up**: nalgebra sparse storage and some
preconditioner internals still retain local `DVector`/`CsrMatrix` conversion
bridges; cfd-validation numerical result/error storage still uses `DVector`;
broad cfd-validation nextest still fails in the existing venturi
cross-fidelity convergence tests.

---

# Sprint 1.96.128 Checklist: Preconditioner Trait Leto Vectors
**Goal**: Move the public `Preconditioner::apply_to` residual/result boundary
from nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `Preconditioner::apply_to` accepts Leto residual and output arrays.
- [x] cfd-math concrete preconditioners implement or call the Leto public
  boundary.
- [x] cfd-math preconditioner tests exercise the Leto boundary.
- [x] cfd-1d network `DiagJacobi` implements the Leto preconditioner contract.
- [x] Targeted scans find no old public `apply_to` DVector signature.

### Closure
- [x] Run cfd-math all-target check/clippy.
- [x] Run cfd-1d lib check/clippy.
- [x] Run focused cfd-math solver/preconditioner nextest (131/131 passed).
- [x] Run cfd-math/cfd-1d fmt check.
- [x] Run targeted `apply_to` DVector-signature residue scan.

**Residual**: `LinearOperator::apply`, nalgebra sparse preconditioner
internals, and iterative solver workspaces still retain `DVector`/`CsrMatrix`
conversion bridges; cfd-validation numerical result/error storage still uses
`DVector`; broad cfd-validation nextest still fails in the existing venturi
cross-fidelity convergence tests.

---

# Sprint 1.96.127 Checklist: Iterative Solver Trait Leto Vectors
**Goal**: Move the public `IterativeLinearSolver::solve` RHS/result boundary
from nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `IterativeLinearSolver::solve` accepts Leto RHS and mutable solution
  arrays.
- [x] CG, BiCGSTAB, and GMRES implement the Leto public boundary.
- [x] cfd-math tests call `solve` through Leto arrays.
- [x] cfd-1d, cfd-2d, and cfd-3d iterative solver callers convert through
  local Leto bridges.
- [x] Targeted scans find no old direct DVector call pattern for `solve`.

### Closure
- [x] Run cfd-math fmt/check/clippy.
- [x] Run focused cfd-math solver nextest (61/61 passed).
- [x] Run cfd-1d/cfd-2d/cfd-3d focused checks and clippy.
- [x] Run cfd-validation no-default all-target clippy.
- [x] Run targeted DVector call-site residue scan.
- [x] Run `git diff --check`.

**Residual**: `LinearOperator::apply`, `Preconditioner::apply_to`,
preconditioners, and internal iterative workspaces still expose nalgebra
`DVector`/`CsrMatrix`; cfd-validation numerical result/error storage still
uses `DVector`; broad cfd-validation nextest still fails in the existing
venturi cross-fidelity convergence tests.

---

# Sprint 1.96.126 Checklist: cfd-math LinearSolver Trait Leto Vectors
**Goal**: Move the public iterative solver `solve_system` vector boundary from
nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `LinearSolver::solve_system` accepts Leto RHS and optional initial-guess
  arrays and returns a Leto result array.
- [x] CG, BiCGSTAB, and GMRES implement the Leto public boundary.
- [x] cfd-validation numerical solver validation calls the Leto public API.
- [x] Targeted scans find no old `solve_system` DVector signature residue.

### Closure
- [x] Run `cargo fmt -p cfd-math -p cfd-validation --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo check -p cfd-validation --no-default-features --lib`.
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `cargo clippy -p cfd-validation --no-default-features --lib -- -D
  warnings`.
- [x] Run `cargo clippy -p cfd-validation --no-default-features --all-targets
  -- -D warnings`.
- [x] Run `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features
  conjugate_gradient bicgstab gmres --status-level fail` (58/58 passed).
- [x] Run targeted `solve_system` signature and DVector-residue scans.
- [x] Run `git diff --check`.

**Residual**: `IterativeLinearSolver`, `LinearOperator`, `Preconditioner`,
preconditioners, and internal iterative workspaces still expose nalgebra
`DVector`/`CsrMatrix`; cfd-validation numerical solver result/error storage
still uses `DVector`; broad cfd-validation nextest still fails in the existing
venturi cross-fidelity convergence tests.

---

# Sprint 1.96.125 Checklist: cfd-validation Leto SpMV and Scalar Bounds
**Goal**: Clear the cfd-validation blockers exposed by the migrated public
Leto SpMV and cfd-1d/cfd-math scalar contracts.

**Success Criteria**:
- [x] cfd-validation SpMV benchmark callers use `leto::Array1`, not nalgebra
  `DVector`, at public `cfd_math::sparse::spmv` call sites.
- [x] linear-solver validation uses `ValidationScalar` for sparse
  linear-operator calls.
- [x] 1D blood-flow literature validations use `ValidationScalar` for cfd-1d
  network solver dispatch.
- [x] `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings` no longer fails in cfd-validation.

### Closure
- [x] Run `cargo fmt -p cfd-validation --check`.
- [x] Run `cargo check -p cfd-validation --no-default-features --lib`.
- [x] Run `cargo clippy -p cfd-validation --no-default-features --lib -- -D
  warnings`.
- [x] Run `cargo clippy -p cfd-validation --no-default-features --all-targets
  -- -D warnings`.
- [x] Run `cargo clippy -p cfd-2d --no-default-features --all-targets -- -D
  warnings`.
- [x] Run cfd-validation targeted SpMV DVector-residue scan.
- [x] Run `cargo nextest run -p cfd-validation --no-default-features
  benchmark --status-level fail` (40/40 passed).
- [x] Run `git diff --check`.

**Residual**: broad `cargo nextest run -p cfd-validation --no-default-features
--status-level fail` fails in two venturi cross-fidelity convergence tests:
`microventuri_35um_case_produces_converged_informative_2d_result` and
`option2_selected_45um_geometry_routes_to_fallback_and_converges`.

---

# Sprint 1.96.124 Checklist: Solver Chain and FEM Leto Vectors
**Goal**: Move the next solver-chain and FEM/direct-fallback vector boundary
from nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `LinearSolverChain::solve` accepts and returns `leto::Array1<T>`.
- [x] `LinearSolverChain::solve_with_guess` accepts Leto RHS and initial guess
  arrays and returns a Leto result array.
- [x] cfd-2d momentum and pressure direct fallbacks route through one private
  Leto bridge into `DirectSparseSolver`.
- [x] cfd-3d FEM and projection sparse assembly route `build_with_rhs` through
  one private FEM Leto bridge.
- [x] cfd-2d/cfd-3d scalar seams carry the Leto real-scalar provider bound.

### Closure
- [x] Run `cargo fmt -p cfd-math -p cfd-1d -p cfd-2d -p cfd-3d --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo check -p cfd-1d --no-default-features --lib`.
- [x] Run `cargo check -p cfd-2d --no-default-features --lib`.
- [x] Run `cargo check -p cfd-3d --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features chain
  direct_solver core_solver simple_gmres --status-level fail` (4/4 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `cargo clippy -p cfd-2d --no-default-features --lib -- -D
  warnings`.
- [x] Run `cargo clippy -p cfd-3d --no-default-features --lib -- -D
  warnings`.
- [x] Run provider-residue scans for direct solver, chain, and `build_with_rhs`
  call sites.

**Residual**: `cargo clippy -p cfd-2d --no-default-features --all-targets --
-D warnings` currently fails in `cfd-validation` because validation still calls
public Leto SpMV with nalgebra vectors and needs Leto scalar bounds propagated
through generic validation/1D-reference paths.

---

# Sprint 1.96.123 Checklist: cfd-math Direct Solver Leto Vectors
**Goal**: Move `DirectSparseSolver` vector inputs and outputs from nalgebra
`DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `DirectSparseSolver::solve` accepts `leto::Array1<T>`.
- [x] `DirectSparseSolver::solve` returns `leto::Array1<T>`.
- [x] rsparse RHS conversion and solution construction operate on Leto arrays.
- [x] Dense fallback uses the existing Leto-array path without a DVector
  wrapper.
- [x] `LinearSolverChain` owns the only remaining direct-tier conversion at
  its current nalgebra boundary.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features direct_solver
  chain core_solver simple_gmres --status-level fail` (4/4 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `cargo doc -p cfd-math --no-default-features --no-deps`.
- [x] Run `cargo test --doc -p cfd-math --no-default-features` (3 passed, 3
  ignored).
- [x] Run targeted direct-solver DVector-signature residue scan.
- [x] Run `git diff --check`.

---

# Sprint 1.96.122 Checklist: cfd-math Sparse Builder Leto RHS
**Goal**: Move sparse matrix Dirichlet RHS assembly from nalgebra `DVector` to
Leto `Array1`.

**Success Criteria**:
- [x] `SparseMatrixBuilder::build_with_rhs` accepts `leto::Array1<T>`.
- [x] Dirichlet column elimination mutates the Leto RHS directly.
- [x] `SparseMatrixBuilder::build` no longer carries a dummy nalgebra RHS.
- [x] Direct-solver and block-preconditioner builder call sites use Leto RHS
  values for assembly.
- [x] Sparse tests assert value-semantic Dirichlet column elimination into the
  Leto RHS.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse
  direct_solver block_preconditioner --status-level fail` (25/25 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run targeted `build_with_rhs` residue scan.
- [x] Run `git diff --check`.

---

# Sprint 1.96.121 Checklist: cfd-math Public SpMV Leto Vectors
**Goal**: Move public sparse SpMV vector arguments from nalgebra `DVector` to
Leto `Array1`.

**Success Criteria**:
- [x] The then-public SpMV entry points accepted `leto::Array1` input/output
  vectors; the redundant parallel-named wrapper was removed on 2026-07-10.
- [x] The nalgebra `DVector` SpMV path is private to `LinearOperator for
  CsrMatrix`.
- [x] Sparse tests, GMRES/AMG integration tests, interpolation quality checks,
  and the SpMV benchmark call public SpMV with Leto arrays.
- [x] Public SpMV signature scan shows no public `DVector` SpMV wrapper.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse spmv
  interpolation amg simple_gmres core_solver --status-level fail` (40/40
  passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `git diff --check` for the touched sparse/solver test/bench files.

---

# Sprint 1.96.120 Checklist: cfd-math SparseMatrixExt Leto Vectors
**Goal**: Move sparse extension diagonal and scaling vector surfaces from
nalgebra `DVector` to Leto `Array1`.

**Success Criteria**:
- [x] `SparseMatrixExt::diagonal` returns `leto::Array1<T>`.
- [x] `SparseMatrixExt::set_diagonal`, `scale_rows`, and `scale_columns`
  accept `leto::Array1<T>`.
- [x] The sparse extension scaling/condition tests use Leto vectors for row
  and column scaling.
- [x] `JacobiPreconditioner::new` consumes the Leto diagonal without restoring
  a nalgebra compatibility vector.
- [x] Targeted residue scan finds no `DVector` signatures for
  `SparseMatrixExt` diagonal/scaling methods.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse basic
  --status-level fail` (21/21 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `git diff --check` for the touched sparse/basic files.

---

# Sprint 1.96.119 Checklist: cfd-math Multigrid Smoother/Cycle Leto Vectors
**Goal**: Move the multigrid smoother and cycle vector paths from nalgebra
`DVector` contracts to Leto `Array1` while keeping the remaining public
sparse/linear-solver boundary explicit.

**Success Criteria**:
- [x] `MultigridLevel`, `AMGHierarchy`, and `MultigridSmoother` use
  `MultigridVector<T> = leto::Array1<T>`.
- [x] Jacobi, Gauss-Seidel, symmetric Gauss-Seidel, SSOR, and Chebyshev
  smoothers operate on Leto vectors.
- [x] Smoother and cycle residuals use `sparse::spmv_array` backed by
  `leto_ops::spmv_into`.
- [x] V-cycle, W-cycle, and F-cycle functions consume and return
  `MultigridVector<f64>`.
- [x] Coarsest cycle solves consume the Leto dense-bridge array path.
- [x] AMG V-cycle recursion uses Leto vectors internally and bridges to
  nalgebra `DVector` only at `Preconditioner::apply_to`.
- [x] `smoothers.rs` and `cycles.rs` have no `DMatrix`, `DVector`,
  `nalgebra::`, `num_traits`, `ndarray`, `rayon`, `tokio`, `rustfft`, or
  `wgpu` residue.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features
  multigrid::cycles smoothers --status-level fail` (10/10 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `git diff --check` for the migrated sparse/multigrid files.

---

# Sprint 1.96.118 Checklist: cfd-math GMG Leto/Eunomia Migration
**Goal**: Move the geometric multigrid cone from nalgebra dense storage and
`num_traits` conversions to Leto arrays and Eunomia scalar traits.

**Success Criteria**:
- [x] `GeometricMultigrid` stores level matrices as `leto::Array2`.
- [x] GMG linear solve, FAS solve, and `NonlinearOperator` surfaces use
  `leto::Array1`.
- [x] Jacobi relaxation, restriction, prolongation, residual, and residual norm
  operate on Leto arrays.
- [x] GMG scalar constants use Eunomia `FloatElement`/`NumericElement`.
- [x] GMG tests use Leto arrays and retain value-semantic matrix, restriction,
  solve, and FAS assertions.
- [x] The GMG directory has no nalgebra/ndarray/num-traits/rayon/tokio/
  rustfft/wgpu residue.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features gmg
  --status-level fail` (5/5 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `git diff --check` for the GMG files.

---

# Sprint 1.96.117 Checklist: cfd-math GMRES Leto/Eunomia Workspace
**Goal**: Move GMRES internal dense Krylov and Givens state from nalgebra
storage to Leto arrays and Eunomia scalar operations.

**Success Criteria**:
- [x] `gmres::arnoldi` stores the Krylov basis and Hessenberg matrix as
  `leto::Array2`.
- [x] `gmres::solver::GMRESWorkspace` stores Krylov basis, Hessenberg matrix,
  Givens coefficients, and least-squares RHS in `leto::Array2`/`Array1`.
- [x] `gmres::givens` uses `leto::Array1`/`Array2` and Eunomia
  `RealField`/`NumericElement`.
- [x] `LinearSolverChain` propagates the Eunomia real-field bound required by
  the migrated GMRES constructor.
- [x] The GMRES module has no `DMatrix`, `ndarray`, `num_traits`,
  `num_complex`, `rayon`, `tokio`, `rustfft`, or `wgpu` residue.
- [x] Remaining nalgebra `DVector`/`RealField` public linear-solver boundary is
  recorded as residual migration work.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features gmres
  --status-level fail` (21/21 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `git diff --check` for the GMRES and chain files.

---

# Sprint 1.96.116 Checklist: cfd-math AMG Restriction Leto Arrays
**Goal**: Move the standalone AMG restriction dense-transfer utilities from
nalgebra dense matrices/vectors to Leto arrays and Leto-ops dense products.

**Success Criteria**:
- [x] Restriction construction functions return `leto::Array2<f64>`.
- [x] Restriction validation reads Leto matrices directly and no longer builds
  nalgebra transpose/difference temporaries.
- [x] `restrict_vector` consumes/returns `leto::Array1<f64>` and asserts
  concrete `P^T v` values in tests.
- [x] `restrict_matrix` delegates `R * A * P` to `leto_ops::MatrixProduct`.
- [x] The migrated file has no nalgebra/ndarray/num-traits/rayon/tokio/
  rustfft/wgpu residue.
- [x] Remaining public sparse/linear-solver nalgebra boundaries are recorded
  as residual migration work.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features restriction
  --status-level fail` (7/7 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.

---

# Sprint 1.96.115 Checklist: cfd-math Linear Solver Leto Dense Bridge
**Goal**: Move legacy CSR-to-dense linear-solver computations from local
nalgebra/manual LU paths to one Leto dense linear-algebra bridge.

**Success Criteria**:
- [x] `linear_solver::dense_bridge` builds `leto::Array2`/`Array1` values from
  the current CSR/vector boundary and solves through `leto_ops::solve`.
- [x] `DirectSparseSolver` dense fallback consumes the bridge.
- [x] Multigrid cycle coarsest small-system solve consumes the bridge.
- [x] The old multigrid local Gaussian-elimination helper is removed.
- [x] The returned solution crosses back to nalgebra `DVector` only at the
  existing public API boundary.
- [x] `LinearSolverChain` carries the Leto real-scalar bound required by the
  direct-solver provider path.
- [x] Remaining public `nalgebra_sparse::CsrMatrix`/`DVector` and `rsparse`
  primary sparse-LU boundaries are recorded as residual migration work.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features
  multigrid::cycles direct_solver --status-level fail` (9/9 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.

---

# Sprint 1.96.114 Checklist: cfd-math Sparse Builder Leto Construction Bridge
**Goal**: Centralize sparse CSR bridge logic and route cfd-math sparse builder
construction through Leto provider validation.

**Success Criteria**:
- [x] Leto/nalgebra CSR conversion lives in one `sparse::bridge` module.
- [x] DVector-to-Leto vector view creation lives in `sparse::bridge`.
- [x] `SparseMatrixBuilder::{build,build_with_rhs,build_parallel}` constructs
  and validates `leto_ops::CsrMatrix` before converting at the legacy boundary.
- [x] `ParallelAssembly::block_diagonal` routes its empty matrix case through
  Leto CSR construction.
- [x] Remaining public `nalgebra_sparse::CsrMatrix`/`DVector` storage boundary
  is recorded as the next larger solver API migration.

### Closure
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse
  --status-level fail` (18/18 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.

---

# Sprint 1.96.113 Checklist: cfd-math Sparse Extension Leto Provider Consumption
**Goal**: Consume Leto-owned CSR utility operations from
`cfd-math::sparse::SparseMatrixExt`.

**Success Criteria**:
- [x] `diagonal`, `scale`, `scale_rows`, `scale_columns`,
  `frobenius_norm`, `is_diagonally_dominant`, and `condition_estimate`
  delegate through `leto_ops::CsrMatrix`.
- [x] `JacobiPreconditioner::new` carries the Leto scalar bound required by
  provider-backed diagonal extraction.
- [x] Sparse tests cover scaling values, shape-error paths, rectangular
  condition-estimate rejection, and the condition-estimate value oracle.
- [x] Remaining sparse/linear-solver `nalgebra_sparse::CsrMatrix`/`DVector`
  storage boundary is recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse
  --status-level fail` (18/18 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.

---

# Sprint 1.96.112 Checklist: cfd-math AMG Leto SpGEMM Consumption
**Goal**: Consume the Leto CSR×CSR product provider in the cfd-math sparse and
AMG Galerkin path.

**Success Criteria**:
- [x] `cfd-math::sparse::try_sparse_sparse_mul` delegates CSR×CSR products to
  `leto_ops::spgemm`.
- [x] `cfd-math::sparse::try_spmv` delegates matrix-vector products to
  `leto_ops::spmv_into`.
- [x] AMG recompute and hierarchy setup use the fallible Leto-backed sparse
  product path.
- [x] AMG restriction construction uses Leto-backed CSR transpose instead of
  `nalgebra_sparse::transpose_as_csc`.
- [x] Sparse product tests cover values, sorted output rows, and shape
  mismatch; sparse transpose tests cover structure and round-trip semantics.
- [x] Remaining sparse/linear-solver `nalgebra_sparse::CsrMatrix`/`DVector`
  storage boundary is recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features sparse
  --status-level fail` (17/17 passed).
- [x] Run `cargo nextest run -p cfd-math --no-default-features interpolation
  --status-level fail` (15/15 passed).
- [x] Run `cargo nextest run -p cfd-math --no-default-features amg
  --status-level fail` (6/6 passed).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.

---

# Sprint 1.96.111 Checklist: Leto-ops CSR Product Provider Gap
**Goal**: Close the upstream CSR×CSR Leto provider gap blocking the
`cfd-math::sparse`/AMG migration.

**Success Criteria**:
- [x] `leto_ops::spgemm` computes CSR×CSR products without a CFDrs-local
  sparse multiply.
- [x] CSR output rows are sorted and exact-zero cancellations are omitted.
- [x] `CsrRow::nnz` exposes row cardinality for sparse-pattern consumers.
- [x] Leto-ops focused sparse tests cover product values, cancellation, and
  shape mismatch.
- [x] Remaining CFDrs sparse/AMG consumer migration is recorded.

### Closure
- [x] Run `cargo fmt -p leto-ops --check`.
- [x] Run `cargo check -p leto-ops`.
- [x] Run `cargo nextest run -p leto-ops --test ops_tests sparse
  --status-level fail` (14/14 passed).
- [x] Run `cargo clippy -p leto-ops --all-targets -- -D warnings`.
- [x] Run `cargo doc -p leto-ops --no-deps`.

---

# Sprint 1.96.110 Checklist: cfd-math SIMD Leto/Eunomia Providers
**Goal**: Remove nalgebra vector/scalar ownership from the cfd-math SIMD cone.

**Success Criteria**:
- [x] `SimdVectorOps` uses Leto `Array1` instead of nalgebra `DVector`.
- [x] SIMD sparse matvec routes through `leto_ops::CsrMatrix`/`spmv`.
- [x] Generic SIMD, field, and vectorization bounds use Eunomia scalar traits.
- [x] SIMD integration tests use Leto arrays and Leto-ops CSR matrices.
- [x] `cargo check -p cfd-math --no-default-features --lib` passes.
- [x] `cargo nextest run -p cfd-math --no-default-features simd` passes.
- [x] Remaining cfd-math sparse/linear-solver provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features simd
  --status-level fail` (26/26 passed, 318 skipped).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run focused SIMD direct-provider residue scan.

---

# Sprint 1.96.109 Checklist: cfd-math Nonlinear Solver Leto Vectors
**Goal**: Remove nalgebra dense vector/matrix ownership from the cfd-math
nonlinear-solver cone.

**Success Criteria**:
- [x] Anderson/JFNK vector math uses Leto `Array1`/`Array2`.
- [x] Shared nonlinear-solver vector helpers own Leto slice, scaling, dot, and
  norm operations.
- [x] JFNK scalar math routes through Eunomia `RealField`/`FloatElement`.
- [x] Nonlinear-solver tests use Leto arrays instead of nalgebra matrices.
- [x] `cargo check -p cfd-math --no-default-features --lib` passes.
- [x] `cargo nextest run -p cfd-math --no-default-features nonlinear`
  passes.
- [x] Remaining cfd-math sparse/linear-solver provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features nonlinear
  --status-level fail` (9/9 passed, 335 skipped).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run focused nonlinear-solver direct-provider residue scan.

---

# Sprint 1.96.108 Checklist: cfd-math DG Leto Dense Arrays
**Goal**: Remove nalgebra dense vector/matrix ownership from the cfd-math
high-order DG cone.

**Success Criteria**:
- [x] `DGBasis` and `DGSolution` use Leto dense arrays.
- [x] Numerical flux, DG operator, limiter, solver, and time-integrator APIs
  use Leto `Array1`/`Array2`.
- [x] DG mass-matrix projection, derivative, RHS, and implicit Newton
  correction solves surface typed Leto solver errors.
- [x] DG docs and DG-related benchmarks compile against Leto arrays.
- [x] `cargo check -p cfd-math --no-default-features --lib` passes.
- [x] DG benchmark target checks pass.
- [x] `cargo nextest run -p cfd-math --no-default-features dg` passes.
- [x] Remaining cfd-math sparse/linear-solver provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo check -p cfd-math --no-default-features --bench
  dg_benchmarks`.
- [x] Run `cargo check -p cfd-math --no-default-features --bench
  flux_alloc_bench`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features dg
  --status-level fail` (62/62 passed, 282 skipped).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run `cargo test --doc -p cfd-math --no-default-features` (3 passed, 3
  ignored).
- [x] Run focused high-order/DG-bench direct-provider residue scan.

---

# Sprint 1.96.107 Checklist: cfd-math Spectral Leto Dense Arrays
**Goal**: Remove nalgebra dense vector/matrix ownership from the cfd-math
high-order spectral cone.

**Success Criteria**:
- [x] `SpectralElement` and `SpectralMesh1D` use Leto dense arrays.
- [x] Spectral differential, interpolation, quadrature, filter, and
  time-integration APIs use Leto `Array1`/`Array2`.
- [x] Spectral assembly local dense matrices/RHS values and dense CSR
  materialization use Leto `Array1`/`Array2`.
- [x] Spectral matrix-vector, stiffness, and dot-product operations share
  local Leto helpers.
- [x] L2 projection solve failures are surfaced as typed errors.
- [x] `cargo check -p cfd-math --no-default-features --lib` passes.
- [x] `cargo nextest run -p cfd-math --no-default-features spectral` passes.
- [x] Remaining cfd-math high-order/sparse/solver provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features spectral
  --status-level fail` (13/13 passed, 331 skipped).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run focused spectral direct-provider residue scan.

---

# Sprint 1.96.106 Checklist: cfd-math WENO Eunomia Scalars
**Goal**: Remove nalgebra scalar trait ownership from the cfd-math high-order
WENO reconstruction cone.

**Success Criteria**:
- [x] `WENO5<T>` and `WENO7<T>` bind on Eunomia scalar provider traits.
- [x] WENO constants use provider-owned conversion helpers.
- [x] WENO nonlinear weights avoid nalgebra math dispatch.
- [x] `cargo check -p cfd-math --no-default-features --lib` passes.
- [x] `cargo nextest run -p cfd-math --no-default-features weno` passes.
- [x] Remaining cfd-math high-order/sparse/solver provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-math --check`.
- [x] Run `cargo check -p cfd-math --no-default-features --lib`.
- [x] Run `cargo nextest run -p cfd-math --no-default-features weno
  --status-level fail` (6/6 passed, 338 skipped).
- [x] Run `cargo clippy -p cfd-math --no-default-features --all-targets -- -D
  warnings`.
- [x] Run focused WENO direct-provider residue scan.

---

# Sprint 1.96.105 Checklist: cfd-1d Solver-Core Eunomia Scalars
**Goal**: Remove direct num-traits scalar bounds and math dispatch from the
cfd-1d primary network solver-core boundary.

**Success Criteria**:
- [x] `NetworkSolveScalar` no longer names `FromPrimitive`, `ToPrimitive`, or
  `num_traits::Float`.
- [x] Anderson acceleration, convergence checks, solver detection, and
  linear-system helpers use Eunomia `FloatElement`/`NumericElement`.
- [x] Solver-core f64 diagnostics use Eunomia scalar extraction.
- [x] `cargo check -p cfd-1d` passes.
- [x] `cargo nextest run -p cfd-1d` passes.
- [x] Remaining vascular/tests/storage provider seams are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-1d --check`.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).
- [x] Run focused solver-core direct-provider residue scan.

---

# Sprint 1.96.104 Checklist: cfd-1d Network Wrapper Eunomia Scalars
**Goal**: Remove direct num-traits scalar construction and nalgebra scalar
bridge conversions from the cfd-1d network wrapper seam.

**Success Criteria**:
- [x] `EdgeProperties::from` uses `SafeFromF64`.
- [x] Network characteristic length uses `SafeFromUsize`.
- [x] Resistance update, hematocrit propagation, validation, and parallel edge
  conductance use Eunomia `NumericElement` for generic absolute values and
  finite checks.
- [x] Blood f64 bridge values use Eunomia `NumericElement::to_f64` rather than
  `nalgebra::try_convert`.
- [x] Adjacent solver bounds compile without reintroducing direct
  `FromPrimitive` into `MatrixAssembler`.
- [x] `cargo check -p cfd-1d` passes.
- [x] `cargo nextest run -p cfd-1d` passes.
- [x] Residual `NetworkSolveScalar` compatibility bounds are recorded.

### Closure
- [x] Run `cargo fmt -p cfd-1d --check`.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).
- [x] Run focused wrapper/matrix/problem direct-provider residue scan.

---

# Sprint 1.96.103 Checklist: cfd-1d Network Blueprint/Sink Eunomia Scalars
**Goal**: Remove direct scalar construction from the canonical cfd-1d
blueprint-to-network entry point while preserving the current wrapper boundary
for the next seam.

**Success Criteria**:
- [x] `network_from_blueprint` constructs blueprint scalar constants through
  `SafeFromF64`/`SafeFromUsize`.
- [x] Generic quadratic-coefficient policy checks use Eunomia
  `NumericElement::abs`.
- [x] `NetworkBuilderSink` propagates the Atlas conversion bounds required by
  the builder.
- [x] `cargo check -p cfd-1d` passes.
- [x] `cargo nextest run -p cfd-1d` passes.
- [x] Residual `domain/network/wrapper.rs` provider work is recorded.

### Closure
- [x] Run `cargo fmt -p cfd-1d`.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).
- [x] Record `Network::new`/`wrapper.rs` as the next network provider seam.

---

# Sprint 1.96.102 Checklist: cfd-1d Domain Components Eunomia Scalars
**Goal**: Remove direct num-traits scalar construction and math dispatch from
the cfd-1d domain-components provider seam.

**Success Criteria**:
- ✅ `Component<T>` pressure-drop math uses Eunomia `NumericElement`.
- ✅ Component factory defaults use the Atlas scalar conversion seam.
- ✅ Channel, membrane, mixer, pump, valve, and sensor implementations no
  longer import direct `num_traits`.
- ✅ `cargo check -p cfd-1d` passes.
- ✅ `cargo nextest run -p cfd-1d` passes.
- ✅ Residual direct provider and clippy debt outside components is recorded.

### Closure
- [x] Run `cargo fmt -p cfd-1d --check`.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).
- [x] Run focused domain-components provider-residue scan.
- [x] Re-run all-target clippy and confirm no touched `channels.rs` lint
  remains.

---

# Sprint 1.96.101 Checklist: cfd-1d Channel/Branching/Analysis Eunomia Scalars
**Goal**: Remove direct num-traits scalar construction and math dispatch from
the next coherent cfd-1d channel, branching, and network-analysis seam.

**Success Criteria**:
- ✅ Channel flow-regime classification and channel solver formulas use
  `SafeFromF64`, `FloatElement`, and `NumericElement`.
- ✅ Branching network solver bounds no longer carry obsolete `FromPrimitive`
  or `ToPrimitive` requirements.
- ✅ Solver-analysis flow, pressure, resistance, performance, and blood-safety
  paths no longer import direct `num_traits`.
- ✅ `cargo check -p cfd-1d` passes.
- ✅ `cargo nextest run -p cfd-1d` passes.
- ✅ Residual direct `num-traits` package areas are recorded separately.

### Closure
- [x] Run touched-file rustfmt through `cargo fmt -p cfd-1d`.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).
- [x] Run focused channel/branching/analyzer and solver-analysis residue scans.
- [x] Record the remaining cfd-1d provider residuals.

---

# Sprint 1.96.100 Checklist: cfd-core Fluid Dynamics Operations Eunomia Scalars
**Goal**: Remove direct num-traits scalar construction from core
fluid-dynamics flow-field operations.

**Success Criteria**:
- ✅ `operations.rs` imports Eunomia `NumericElement` instead of direct
  `num_traits::FromPrimitive`.
- ✅ Vorticity and divergence finite-difference constants use provider-owned
  zero/one/two identities.
- ✅ Kinetic-energy and enstrophy half factors use provider-owned scalar
  construction.
- ✅ Focused operations value tests pass under nextest.
- ✅ The remaining nalgebra vector/storage migration is recorded separately.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core fluid_dynamics::operations`.
- [x] Run focused `operations.rs` scalar-provider residue scan.
- [x] Record the remaining fluid-dynamics nalgebra storage residuals.

---

# Sprint 1.96.99 Checklist: cfd-core Fluid Dynamics Service Eunomia Scalars
**Goal**: Remove direct num-traits scalar construction and math dispatch from
core fluid-dynamics pipe-flow service formulas.

**Success Criteria**:
- ✅ `service.rs` imports Eunomia `FloatElement`/`NumericElement` instead of
  direct `num_traits::{Float, FromPrimitive}`.
- ✅ Pipe pressure-drop and friction-factor constants use Eunomia scalar
  construction without silent fallback conversion.
- ✅ Colebrook-White and Haaland math uses Eunomia powers, square roots,
  logarithms, and absolute convergence checks.
- ✅ Focused service tests pass under nextest.
- ✅ The remaining `operations.rs` and nalgebra storage residuals are recorded
  separately.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core fluid_dynamics::service`.
- [x] Run focused `service.rs` scalar-provider residue scan.
- [x] Record the remaining fluid-dynamics provider residuals.

---

# Sprint 1.96.98 Checklist: cfd-core Flow Regime Eunomia Scalars
**Goal**: Remove direct nalgebra/num-traits scalar conversion from core
flow-regime classification.

**Success Criteria**:
- ✅ `flow_regimes.rs` imports Eunomia `RealField`/`NumericElement` instead of
  direct `nalgebra::RealField` or `num_traits::ToPrimitive`.
- ✅ Reynolds, Mach, and combined classification use
  `NumericElement::to_f64` without silent `unwrap_or(0.0)` conversion
  fallback.
- ✅ `FluidDynamicsService::flow_regime` exposes the migrated Eunomia scalar
  contract.
- ✅ Value-semantic threshold tests pass under nextest.
- ✅ The broader service friction-factor and vector-storage residuals are
  recorded separately.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core flow_regime`.
- [x] Run focused `flow_regimes.rs` scalar-provider residue scan.
- [x] Record the remaining fluid-dynamics provider residuals.

---

# Sprint 1.96.97 Checklist: cfd-3d Spectral Diagnostics Eunomia Conversion
**Goal**: Remove direct num-traits scalar conversion from the Apollo/Leto
3D spectral diagnostics path.

**Success Criteria**:
- ✅ `spectral::diagnostics` imports Eunomia `NumericElement` instead of
  direct `num_traits` conversion traits.
- ✅ Kinetic-energy and enstrophy spectra continue to stage Leto arrays for
  Apollo FFTs.
- ✅ Existing diagnostics value tests pass under nextest.
- ✅ Focused residue scan is clean for direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `.to_f64()`, and obsolete conversion-context plumbing.
- ✅ The remaining `VelocityField<T>` nalgebra storage boundary is recorded as
  residual upstream work.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run `cargo nextest run -p cfd-3d diagnostics`.
- [x] Run focused `diagnostics.rs` scalar-provider residue scan.
- [x] Record the broader `VelocityField<T>` and spectral dense-linalg
  residuals.

---

# Sprint 1.96.96 Checklist: cfd-2d MRT/Carreau-Yasuda Eunomia Scalars
**Goal**: Close the remaining LBM MRT and Carreau-Yasuda direct
nalgebra/num-traits scalar-provider residue.

**Success Criteria**:
- ✅ MRT relaxation matrices, moment transforms, and equilibrium moments use
  Eunomia-backed scalar construction and identities.
- ✅ The local Carreau-Yasuda rheology model and BGK collision operator use
  Eunomia `FloatElement`/`NumericElement` APIs.
- ✅ Focused LBM and Carreau-Yasuda value tests pass under nextest.
- ✅ Focused residue scans are clean for direct `RealField`, `num_traits`,
  `FromPrimitive`, `T::from_f64`, `T::zero`, `T::one`, `from_i32`, and
  `Float::` tokens in the migrated files.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo nextest run -p cfd-2d lbm`.
- [x] Run `cargo nextest run -p cfd-2d carreau_yasuda`.
- [x] Run focused MRT/Carreau-Yasuda scalar-provider residue scan.

---

# Sprint 1.96.95 Checklist: cfd-2d LBM Eunomia Scalar Seam
**Goal**: Remove direct nalgebra/num-traits scalar construction from the
authoritative 2D LBM macroscopic, equilibrium, collision trait, and BGK seam.

**Success Criteria**:
- ✅ `macroscopic.rs`, `lattice.rs`, `collision/traits.rs`, and
  `collision/bgk.rs` no longer import or bound on direct `nalgebra`,
  `RealField`, or `num_traits`.
- ✅ LBM constants, zero/one identities, weights, and BGK viscosity/omega
  construction route through Eunomia-backed scalar helpers.
- ✅ Existing focused LBM value tests pass under nextest.
- ✅ Residual MRT/Carreau-Yasuda inherent legacy scalar work is recorded as
  out of scope for this slice and closed by Sprint 1.96.96.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo nextest run -p cfd-2d lbm`.
- [x] Run focused LBM scalar-provider residue scan.
- [x] Record the broader MRT/Carreau-Yasuda and CFDrs provider residuals.

---

# Sprint 1.96.94 Checklist: cfd-3d Wall Functions Eunomia Scalars
**Goal**: Remove direct nalgebra/num-traits scalar math from the 3D Spalding
wall-law helper.

**Success Criteria**:
- ✅ `wall_functions.rs` imports Eunomia scalar traits instead of nalgebra or
  num-traits.
- ✅ Spalding constants and transcendental math route through
  `FloatElement`/`NumericElement`.
- ✅ Existing wall-law value tests pass under nextest.
- ✅ Focused residue scan is clean for direct `nalgebra`/`num_traits` imports
  and bounds.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run `cargo nextest run -p cfd-3d wall_functions`.
- [x] Run focused `wall_functions.rs` scalar-provider residue scan.
- [x] Record the direct WGPU version-boundary residual for a later
  Hephaestus context migration.

---

# Sprint 1.96.93 Checklist: cfd-3d Apollo Fourier Eunomia Bounds
**Goal**: Remove direct num-traits scalar conversion from the Apollo-backed
3D Fourier wrapper and repair downstream provider-bound fallout.

**Success Criteria**:
- ✅ `cfd-3d::spectral::fourier` uses Eunomia conversion traits.
- ✅ Apollo FFT + Leto array execution remains the Fourier path.
- ✅ Downstream blood/cavitation consumers compile with migrated
  `FloatElement` contracts.
- ✅ Focused Fourier integration tests pass under nextest.

### Closure
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run `cargo check -p cfd-validation`.
- [x] Run `cargo nextest run -p cfd-3d --test fourier_validation`.
- [x] Run focused `cfd-3d` Fourier/DES scalar-provider residue scan.

---

# Sprint 1.96.92 Checklist: cfd-core Fåhræus-Lindqvist Eunomia Scalars
**Goal**: Remove the final local blood-model direct num-traits scalar holdout.

**Success Criteria**:
- ✅ `FahraeuasLindqvist<T>::new` uses Eunomia scalar construction.
- ✅ Pries/Secomb relative-viscosity formulas use Eunomia `powf`/`exp`.
- ✅ Relative-viscosity clamping and tube hematocrit use Eunomia identities,
  `abs`, and max dispatch.
- ✅ Existing Fåhræus-Lindqvist and blood-dispatch value tests remain green.
- ✅ Focused and broader blood residue scans are clean for generic direct
  num-traits/math residue.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm Fåhræus-Lindqvist needs scalar constants, exponentials, real
  powers, absolute value, and max clamping.
- [x] Confirm the broader `Fluid<T>` trait still owns the inherited
  `RealField` boundary.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` construction with Eunomia
  `FloatElement`.
- [x] Replace direct zero/one, `powf`, `exp`, `abs`, and max calls with
  Eunomia APIs.

### Phase 3: Closure (50%+)
- [x] Run focused `fahraeus_lindqvist.rs` residue scan.
- [x] Run broader blood residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core fahraeus_lindqvist`.
- [x] Run `cargo nextest run -p cfd-core blood`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.91 Checklist: cfd-core Casson/Carreau Blood Eunomia Scalars
**Goal**: Remove Casson/Carreau/BloodModel direct num-traits scalar dispatch.

**Success Criteria**:
- ✅ `CassonBlood<T>` constructors and scalar math use Eunomia.
- ✅ `CarreauYasudaBlood<T>` remains on Eunomia and compiles through the shared
  dispatch enum.
- ✅ `BloodModel<T>` no longer imports or requires `num_traits::FromPrimitive`.
- ✅ Existing Casson, Carreau-Yasuda, and blood-dispatch value tests remain
  green.
- ✅ Focused touched-file residue scan is clean except for the concrete `f64`
  temperature helper.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm Casson needs scalar constants, exponential temperature scaling,
  square roots, and validation identities.
- [x] Confirm the broader `Fluid<T>` trait still owns the inherited
  `RealField` boundary.

### Phase 2: Execution (10-50%)
- [x] Replace Casson direct `FromPrimitive` construction with Eunomia
  `FloatElement`.
- [x] Replace Casson direct zero/one, `sqrt`, and generic `exp` calls with
  Eunomia APIs.
- [x] Replace `BloodModel` dispatch bound with Eunomia `FloatElement`.

### Phase 3: Closure (50%+)
- [x] Run focused touched-file residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core casson`.
- [x] Run `cargo nextest run -p cfd-core carreau_yasuda`.
- [x] Run `cargo nextest run -p cfd-core blood`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.90 Checklist: cfd-core Cross Blood Eunomia Scalars
**Goal**: Remove Cross blood model-local direct num-traits scalar construction.

**Success Criteria**:
- ✅ `CrossBlood<T>::normal_blood` uses Eunomia scalar construction.
- ✅ `CrossBlood<T>::apparent_viscosity` uses Eunomia zero/one and `powf`.
- ✅ Existing Cross and blood-dispatch value tests remain green.
- ✅ Focused `cross.rs` residue scan is clean.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm Cross blood needs scalar constants, comparisons, and real-power
  evaluation.
- [x] Confirm the broader `Fluid<T>` trait still owns the inherited
  `RealField` boundary.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` construction with Eunomia
  `FloatElement`.
- [x] Replace direct zero/one and generic `powf` calls with Eunomia APIs.

### Phase 3: Closure (50%+)
- [x] Run focused `cross.rs` residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core cross`.
- [x] Run `cargo nextest run -p cfd-core blood`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.89 Checklist: cfd-core Cavitation Scalar Closeout
**Goal**: Remove the remaining cavitation-local nalgebra scalar contracts.

**Success Criteria**:
- ✅ `CavitationModel<T>` and `ZgbParams<T>` use Eunomia scalar bounds.
- ✅ The heterogeneous legacy adapter no longer exposes a fake generic
  `RealField` surface.
- ✅ New value-semantic tests cover mass-transfer and adapter behavior.
- ✅ Focused cavitation residue scan is clean.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm mass-transfer models only need scalar arithmetic, comparisons,
  constants, powers, roots, absolute value, and min/max.
- [x] Confirm the heterogeneous legacy helper already computes through the
  concrete `f64` selective-cavitation model.

### Phase 2: Execution (10-50%)
- [x] Replace `CavitationModel<T: RealField + Copy>` and `ZgbParams<T>` with
  Eunomia `FloatElement` bounds.
- [x] Replace mass-transfer constants and math with Eunomia APIs.
- [x] Replace the fake-generic heterogeneous legacy adapter with an honest
  `f64` contract.
- [x] Add value-semantic regression tests for migrated formulas and adapter
  parity.

### Phase 3: Closure (50%+)
- [x] Run focused cavitation residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core cavitation`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.88 Checklist: cfd-core Nuclei Transport Eunomia Scalars
**Goal**: Remove the nuclei transport-local nalgebra scalar contract.

**Success Criteria**:
- ✅ `nuclei_adjusted_vapor_pressure<T>` uses Eunomia scalar bounds.
- ✅ `NucleiTransportConfig<T>` uses Eunomia scalar bounds and defaults.
- ✅ `NucleiTransport<T>` uses Eunomia scalar bounds.
- ✅ Existing nuclei transport value tests remain green.
- ✅ Focused nuclei transport residue scan is clean.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm nuclei transport only needs scalar arithmetic, comparisons,
  constants, and exponential evaluation.
- [x] Confirm existing tests cover affine pressure coupling, dissolution,
  exponential transit decay, and diffusion accessor behavior.

### Phase 2: Execution (10-50%)
- [x] Replace `nuclei_adjusted_vapor_pressure<T: RealField + Copy>` with
  Eunomia `FloatElement`.
- [x] Replace `NucleiTransportConfig<T: RealField + Copy>` with Eunomia
  `FloatElement`.
- [x] Replace `NucleiTransport<T: RealField + Copy>` with Eunomia
  `FloatElement`.
- [x] Replace direct scalar constants and exponentials with Eunomia APIs.

### Phase 3: Closure (50%+)
- [x] Run focused nuclei transport residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core cavitation`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.87 Checklist: cfd-core Venturi Cavitation Eunomia Scalars
**Goal**: Remove the Venturi cavitation-local nalgebra scalar contract and
direct num-traits construction.

**Success Criteria**:
- ✅ `VenturiCavitation<T>` uses Eunomia scalar bounds.
- ✅ Venturi constants and scalar identities route through Eunomia.
- ✅ Existing Venturi value tests remain green.
- ✅ Focused Venturi residue scan is clean for production code.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm Venturi only needs scalar arithmetic, comparisons, integer
  powers, tangent, absolute value, and constants.
- [x] Confirm existing tests cover the migrated physical formulas.

### Phase 2: Execution (10-50%)
- [x] Replace `VenturiCavitation<T: RealField + Copy>` with
  `VenturiCavitation<T: FloatElement + Copy>`.
- [x] Replace direct `T::zero()`/`T::one()` and `FromPrimitive` construction
  with Eunomia `NumericElement`/`FloatElement`.
- [x] Replace direct generic `powi`, `powf`, `tan`, and `abs` calls with
  explicit Eunomia dispatch.

### Phase 3: Closure (50%+)
- [x] Run focused Venturi/migrated-cone residue scans.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core cavitation`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.86 Checklist: cfd-core Cavitation Eunomia Scalar Cone
**Goal**: Remove nalgebra scalar contracts and direct num-traits construction
from the current touched cavitation scalar cone.

**Success Criteria**:
- ✅ Rayleigh-Plesset and sonoluminescence use Eunomia scalar bounds.
- ✅ Cavitation biological damage uses Eunomia scalar bounds.
- ✅ Cavitation regime classification/analysis uses Eunomia scalar bounds.
- ✅ Cavitation number and material-damage models use Eunomia scalar bounds.
- ✅ New closed-form tests cover the newly migrated number/damage APIs.
- ✅ Focused cavitation residue scan is clean for migrated files.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm touched cavitation files only need scalar arithmetic,
  comparisons, finite checks, powers, roots, exponentials, and constants.
- [x] Identify remaining cavitation holdouts before claiming scope closure.

### Phase 2: Execution (10-50%)
- [x] Replace Rayleigh-Plesset, biological-damage, and regime-analysis bounds
  with Eunomia `FloatElement`/`NumericElement`.
- [x] Replace `CavitationNumber<T>` bounds and constants with Eunomia.
- [x] Replace `CavitationDamage<T>` bounds and constants with Eunomia.
- [x] Add value-semantic closed-form tests for cavitation number and damage.

### Phase 3: Closure (50%+)
- [x] Run focused cavitation residue scan over migrated files.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core cavitation`.
- [x] Run touched-file `git diff --check`.
- [x] Record unrelated full-clippy blocker explicitly.

---

# Sprint 1.96.85 Checklist: cfd-core Hemolysis Eunomia Scalars
**Goal**: Remove the hemolysis-local nalgebra scalar contract and direct
num-traits construction.

**Success Criteria**:
- ✅ `HemolysisCalculator<T>` uses Eunomia scalar bounds.
- ✅ `PlateletActivation<T>` uses Eunomia scalar bounds.
- ✅ Hemolysis constants and zero/one identities route through Eunomia.
- ✅ Focused hemolysis residue scan is clean.
- ✅ Focused compile, nextest, formatting, and diff checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm hemolysis only needs scalar construction, comparisons,
  arithmetic, and exponential evaluation.
- [x] Confirm no active scheduler or ndarray dependency slice is present in
  this turn's exact scan.

### Phase 2: Execution (10-50%)
- [x] Replace `HemolysisCalculator<T: RealField + Copy>` with
  `HemolysisCalculator<T: FloatElement + Copy>`.
- [x] Replace `PlateletActivation<T: RealField + Copy>` with
  `PlateletActivation<T: FloatElement + Copy>`.
- [x] Replace direct `T::zero()`/`T::one()` and `FromPrimitive` construction
  with Eunomia `NumericElement`/`FloatElement`.

### Phase 3: Closure (50%+)
- [x] Run focused hemolysis residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core hemolysis`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.84 Checklist: cfd-2d Grid/FVM Leto-Eunomia Boundary
**Goal**: Remove the grid-level nalgebra vector/scalar contract that FVM still
inherited through `StructuredGrid2D<T>`.

**Success Criteria**:
- ✅ `Grid2D::cell_center` returns `leto::geometry::Vector2`.
- ✅ `StructuredGrid2D` no longer requires `nalgebra::RealField`.
- ✅ Structured grid spacing and cell-center scalar construction use Eunomia.
- ✅ `AdaptiveGrid2D::effective_resolution` uses Eunomia index conversion.
- ✅ `FvmSolver` no longer carries an inherited `RealField` bound.
- ✅ Validation/FDM/LBM consumers compile against Leto vector indexing.
- ✅ Focused grid/FVM residue scans are clean.
- ✅ Focused compile and nextest checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm FVM's remaining `RealField` owner was the structured-grid
  contract.
- [x] Identify grid-center consumers that assume nalgebra `.x`/`.y` fields.

### Phase 2: Execution (10-50%)
- [x] Change `Grid2D::cell_center` to Leto `Vector2`.
- [x] Replace `StructuredGrid2D` spacing and center conversions with Eunomia
  helpers.
- [x] Replace unstructured grid center storage with Leto `Vector2`.
- [x] Replace adaptive effective-resolution factor conversion with Eunomia.
- [x] Update FVM and grid-center consumers for Leto vector indexing.

### Phase 3: Closure (50%+)
- [x] Run focused grid/FVM residue scans.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo check -p cfd-validation`.
- [x] Run `cargo nextest run -p cfd-2d --lib grid`.
- [x] Run `cargo nextest run -p cfd-2d --lib fvm`.
- [x] Run touched-file `git diff --check`.

---

# Sprint 1.96.83 Checklist: cfd-2d FVM Leto Vector Boundary
**Goal**: Replace FVM-local nalgebra `Vector2` storage and operations with the
Atlas Leto geometry provider.

**Success Criteria**:
- ✅ Leto exposes `leto::geometry::Vector2`.
- ✅ Leto fixed vectors provide generic norm and normalization operations.
- ✅ `Face<T>` uses `leto::geometry::Vector2<T>` for centers and normals.
- ✅ `FvmSolver::solve` accepts velocity fields backed by Leto vectors.
- ✅ No FVM-local nalgebra `Vector2` import remains.
- ✅ `cargo check -p leto` passes.
- ✅ `cargo check -p cfd-2d` passes.
- ✅ Focused Leto and FVM nextest checks pass.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm FVM needs only fixed 2-D vector construction, indexing, norm,
  normalization, and dot products.
- [x] Keep `StructuredGrid2D<T>` as the explicit residual grid-level
  `RealField` owner.

### Phase 2: Execution (10-50%)
- [x] Add provider-owned `Vector2<T>` alias in Leto.
- [x] Add generic fixed-vector norm and normalization methods in Leto.
- [x] Replace FVM face geometry and velocity fields with Leto vectors.
- [x] Add direct `leto` dependency to `cfd-2d`.

### Phase 3: Closure (50%+)
- [x] Run focused FVM provider-residue scans.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p leto`.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run focused Leto and FVM `cargo nextest` checks.

---

# Sprint 1.96.82 Checklist: cfd-2d FVM Eunomia Scalar Constants
**Goal**: Remove direct num-traits scalar construction from FVM configuration,
face-center construction, and flux constants.

**Success Criteria**:
- ✅ `FvmConfig` no longer imports `nalgebra::RealField`.
- ✅ `FvmConfig::default` no longer requires `num_traits::FromPrimitive`.
- ✅ Default scalar constants use Eunomia `FloatElement`.
- ✅ Default constants have value-semantic unit coverage.
- ✅ `FvmSolver` declares the required `FloatElement` owner bound where it
  stores `FvmConfig<T>`.
- ✅ `FvmSolver` face centers use Eunomia-backed index conversion.
- ✅ FVM flux schemes no longer require `num_traits::FromPrimitive`.
- ✅ Non-finite flux diffusion coefficients have value-semantic error coverage.
- ✅ `cargo check -p cfd-2d` passes.
- ✅ `cargo check -p cfd-3d` passes through the migrated spectral dependency
  chain.
- ✅ `cargo nextest run -p cfd-2d --lib fvm` passes focused FVM coverage.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `FvmConfig` storage/defaults are a local configuration boundary.
- [x] Keep `FvmSolver`/flux `RealField` residuals explicit.

### Phase 2: Execution (10-50%)
- [x] Replace config `RealField`/`FromPrimitive` bounds with Eunomia
  `FloatElement`.
- [x] Replace default scalar construction with `FloatElement::from_f64`.
- [x] Add value-semantic default-constant coverage.
- [x] Propagate the config owner bound to `FvmSolver`.
- [x] Replace solver face-center conversions with Eunomia helpers.
- [x] Replace flux scalar constants with Eunomia helpers.
- [x] Add invalid-diffusion flux tests.

### Phase 3: Closure (50%+)
- [x] Run focused FVM residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run focused `cargo nextest` for the FVM cone.

---

# Sprint 1.96.81 Checklist: cfd-2d CFL Eunomia Scalars
**Goal**: Remove nalgebra and num-traits scalar contracts from the standalone
2D CFL stability calculator.

**Success Criteria**:
- ✅ `CFLCalculator` no longer imports `nalgebra::RealField`.
- ✅ `CFLCalculator` no longer imports `num_traits::FromPrimitive`.
- ✅ CFL constants, absolute values, and scalar construction use Eunomia.
- ✅ Existing CFL tests pass under `cargo nextest run`.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `CFLCalculator` is a local one-file stability utility.
- [x] Verify Eunomia owns the required constants, `abs`, and scalar
  construction APIs.

### Phase 2: Execution (10-50%)
- [x] Replace `RealField`/`FromPrimitive` bounds with Eunomia bounds.
- [x] Replace `T::zero`, `T::one`, `T::from_f64`, and `.abs()` usage with
  Eunomia calls.
- [x] Replace `min` with an explicit comparison under the reduced trait bound.

### Phase 3: Closure (50%+)
- [x] Run focused CFL residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo nextest run -p cfd-2d --lib cfl` (9/9 passed).
- [x] Record unrelated broad-filter integration-test bound debt.

---

# Sprint 1.96.80 Checklist: cfd-core Material Eunomia Traits
**Goal**: Remove nalgebra scalar trait contracts from cfd-core material
solid/interface properties while preserving the fluid database boundary.

**Success Criteria**:
- ✅ `SolidProperties` and `InterfaceProperties` no longer require
  `nalgebra::RealField`.
- ✅ `ElasticSolid`, `WettingProperties`, and `FluidSolidInterface` use
  Eunomia scalar bounds.
- ✅ `MaterialDatabase` keeps `RealField` only through `Fluid<T>` storage.
- ✅ Value-semantic material tests cover the migrated formulas/constants.

### Phase 1: Foundation & Specs (0-10%)
- [x] Separate material solid/interface contracts from the still-fluid-owned
  `MaterialDatabase` boundary.

### Phase 2: Execution (10-50%)
- [x] Replace material solid/interface `RealField` bounds with
  `FloatElement`/`NumericElement`.
- [x] Route shear modulus and adhesion energy through Eunomia operations.
- [x] Add value-semantic material tests.

### Phase 3: Closure (50%+)
- [x] Run focused material residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core material` (4/4 passed).
- [x] Record residual fluid/material database nalgebra coupling.

---

# Sprint 1.96.79 Checklist: cfd-core Velocity Leto Vector
**Goal**: Replace the immediate cfd-core velocity vector value object with the
Leto provider vector while preserving serialization.

**Success Criteria**:
- ✅ `Velocity` stores `leto::geometry::Vector3`.
- ✅ `PhysicalParameters::gravity` stores `leto::geometry::Vector3`.
- ✅ Leto geometry supports Serde for CFDrs serialized value objects.
- ✅ `cfd-core` compiles with direct Leto consumption.

### Phase 1: Foundation & Specs (0-10%)
- [x] Verify Leto exposes the needed 3D vector operations.
- [x] Identify Serde as the provider-side missing contract.

### Phase 2: Execution (10-50%)
- [x] Add Serde derives to Leto fixed 3D geometry value types.
- [x] Add direct Leto dependency to `cfd-core`.
- [x] Switch `Velocity` and `PhysicalParameters::gravity` to Leto `Vector3`.

### Phase 3: Closure (50%+)
- [x] Run focused provider-residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --lib` (183/183 passed).
- [x] Run downstream compile checks for `cfd-2d`, `cfd-3d`, and
  `cfd-validation`.
- [x] Record residual Domain/fluid/material/hemolysis nalgebra work.

---

# Sprint 1.96.78 Checklist: cfd-core Physics Value Eunomia Scalars
**Goal**: Remove direct nalgebra scalar contracts from scalar-only physics
value wrappers while preserving vector-backed velocity for the later Leto/Gaia
slice.

**Success Criteria**:
- ✅ `Temperature`, `Pressure`, `ReynoldsNumber`, and `DimensionlessNumber` no
  longer import `nalgebra::RealField`.
- ✅ Scalar zero, absolute value, and square root use Eunomia
  `NumericElement`.
- ✅ Immediate aggregate owners compile with the propagated Eunomia bounds.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify scalar-only value wrappers versus vector-backed `Velocity`.

### Phase 2: Execution (10-50%)
- [x] Replace scalar wrapper bounds with `FloatElement`/`NumericElement`.
- [x] Propagate bounds through `PhysicalParameters`, `ProblemAggregate`,
  `InitialConditions`, and `SimulationAggregate`.

### Phase 3: Closure (50%+)
- [x] Run focused scalar-provider residue scan.
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --lib` (183/183 passed).
- [x] Record residual vector/material/hemolysis nalgebra migration work.

---

# Sprint 1.96.77 Checklist: cfd-math WGPU Boundary Reduction
**Goal**: Remove the direct `wgpu` dependency from `cfd-math` and route GPU
metric synchronization through `cfd-core` while preserving existing kernels.

**Success Criteria**:
- ✅ `cfd-math` no longer declares optional `wgpu`.
- ✅ `cfd-math/gpu` no longer enables `dep:wgpu`.
- ✅ The root package `gpu` feature no longer enables `dep:wgpu`.
- ✅ `cfd-math::linear_solver::operators::gpu` has no direct `wgpu` import.
- ✅ `GpuContext` exposes synchronization and timestamp-query capability.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm the only `cfd-math` WGPU usage is metric synchronization and
  timestamp-query capability.

### Phase 2: Execution (10-50%)
- [x] Add `GpuContext::synchronize`.
- [x] Add `GpuContext::supports_timestamp_queries`.
- [x] Remove direct `wgpu` dependency and feature activation from `cfd-math`
  and root package features.

### Phase 3: Closure (50%+)
- [x] Run focused WGPU residue scan.
- [x] Run `cargo check -p cfd-math --features gpu`.
- [x] Run `cargo nextest run -p cfd-math --features gpu --lib` (298/298
  passed).
- [x] Run root GPU feature check: `cargo check -p cfd-suite --features gpu`.
- [x] Run touched-file rustfmt and `git diff --check`.
- [x] Record residual raw WGPU kernel/buffer/pipeline ownership in `cfd-core`.

---

# Sprint 1.96.76 Checklist: cfd-core Hephaestus GPU Probe
**Goal**: Route `cfd-core` GPU availability detection through Hephaestus device
acquisition while preserving the existing raw WGPU kernels for the next
provider-migration slice.

**Success Criteria**:
- ✅ `cfd-core` depends on optional `hephaestus-wgpu` under the `gpu` feature.
- ✅ `ComputeBackend::detect_gpu_support` uses
  `hephaestus_wgpu::WgpuDevice::try_default`.
- ✅ `cargo tree --workspace -i hephaestus-wgpu` shows the active
  `hephaestus-wgpu -> cfd-core` provider path.
- ✅ `cargo check -p cfd-core --features gpu` passes.
- ✅ `cargo nextest run -p cfd-core --features gpu --lib` passes.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `rustfft` is absent from the workspace graph.
- [x] Confirm `ndarray` remains only through `numpy -> cfd-python`.
- [x] Confirm the active runtime GPU gap is direct WGPU probing and raw GPU
  ownership in `cfd-core`/`cfd-math`.

### Phase 2: Execution (10-50%)
- [x] Add `hephaestus-wgpu` as the optional cfd-core GPU provider dependency.
- [x] Replace direct raw WGPU adapter probing with Hephaestus device
  acquisition.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --features gpu`.
- [x] Run `cargo nextest run -p cfd-core --features gpu --lib` (183/183
  passed).
- [x] Run `cargo tree --workspace -i hephaestus-wgpu`.
- [x] Record residual raw WGPU buffer/pipeline/shader migration risk.

---

# Sprint 1.96.75 Checklist: cfd-1d Non-Python ndarray Path Removal
**Goal**: Remove the stale `sprs -> ndarray` dependency path from cfd-1d and
drop the unused root workspace `ndarray` declaration.

**Success Criteria**:
- ✅ `cfd-1d` no longer declares `sprs`.
- ✅ The root workspace no longer declares unused `ndarray`.
- ✅ `cargo tree -p cfd-1d -i ndarray` reports no matching package.
- ✅ `cargo tree -p cfd-3d -i ndarray` reports no matching package.
- ✅ cfd-1d tests compile against the migrated Eunomia fluid constructor bounds.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `sprs` is manifest-only residue in cfd-1d.
- [x] Confirm `nalgebra-sparse` is still actively used and out of this slice.

### Phase 2: Execution (10-50%)
- [x] Remove `sprs` from cfd-1d.
- [x] Remove unused workspace `ndarray`.
- [x] Add direct cfd-1d Eunomia dependency and test helper bounds.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo update -p sprs`.
- [x] Run cfd-1d/cfd-3d `ndarray` inverse-tree checks.
- [x] Run `cargo check -p cfd-1d`.
- [x] Run `cargo nextest run -p cfd-1d` (725/725 passed, 3 skipped).

---

# Sprint 1.96.74 Checklist: cfd-3d Apollo ndarray Removal
**Goal**: Resolve CFDrs' Apollo packages through the side-by-side Atlas Apollo
checkout and remove active ndarray-backed Apollo package resolution from the
cfd-3d graph.

**Success Criteria**:
- ✅ `apollo-fft` and `apollo-nufft` resolve from `D:\atlas\repos\apollo`.
- ✅ Active Apollo package trees have no `ndarray` matches.
- ✅ Apollo source/manifest/lock scan has no `ndarray` hits.
- ✅ cfd-3d FEM and spectral call sites compile against Eunomia scalar/complex
  contracts.
- ✅ cfd-validation dev-dependency code carries the required `FloatElement`
  bounds for cfd-3d nextest.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit active Apollo package resolution and remaining `ndarray` inverse
  dependency path.

### Phase 2: Execution (10-50%)
- [x] Patch CFDrs to side-by-side Atlas Apollo and provider checkouts.
- [x] Update cfd-3d FEM/spectral call sites for Eunomia provider contracts.
- [x] Propagate validation `FloatElement` bounds and conversion dispatch.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run `cargo nextest run -p cfd-3d` (394/394 passed; one slow test at
  16.394s).
- [x] Verify active Apollo trees and Apollo source/manifest/lock scan have no
  `ndarray` hits.

---

# Sprint 1.96.73 Checklist: cfd-math Eunomia Geometric Multigrid
**Goal**: Remove direct scalar conversion provider residue from geometric
multigrid and clear the remaining AMG `FromPrimitive` bound.

**Success Criteria**:
- ✅ GMG hierarchy construction uses Eunomia scalar conversion.
- ✅ GMG Poisson matrix constants use Eunomia scalar conversion.
- ✅ GMG transfer weights use Eunomia scalar conversion.
- ✅ AMG no longer carries direct `num_traits::FromPrimitive` bounds/imports.
- ✅ Value-semantic tests cover Poisson stencil constants and full-weighting
  restriction values.
- ✅ `multigrid` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  direct `T::from_f64`, direct `T::from_usize`, conversion fallback,
  `from_f64_or`, `SafeFromF64`, stale `rayon`/`tokio`, `rustfft`, `ndarray`, or
  `num_complex` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit remaining provider residue in `multigrid`.
- [x] Confirm the bounded slice preserves current nalgebra dense/sparse/vector APIs.

### Phase 2: Execution (10-50%)
- [x] Replace GMG hierarchy and Poisson scalar constants with Eunomia.
- [x] Replace GMG transfer weights with Eunomia.
- [x] Remove stale AMG `FromPrimitive` import and bounds.
- [x] Add value-semantic Poisson stencil and restriction-weight tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run multigrid-wide residue scan.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math geometric_multigrid poisson_matrix restrict_residual fas_solve amg` (11/11 passed).

---

# Sprint 1.96.72 Checklist: cfd-math Eunomia Multigrid Interpolation
**Goal**: Remove direct scalar conversion provider residue from multigrid
interpolation operators and interpolation quality metrics.

**Success Criteria**:
- ✅ Interpolation operators use Eunomia `FloatElement` for scalar constants and
  index-distance conversion.
- ✅ Interpolation quality metrics use Eunomia `NumericElement` for row-sum
  extraction, absolute values, and sparsity ratios.
- ✅ Touched tests no longer print debug output or use direct index-distance
  numeric casts.
- ✅ `multigrid/interpolation.rs` has no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, conversion
  fallback, `from_f64_or`, `SafeFromF64`, stale `rayon`, direct `as f64`, or
  `.to_f64()` fallback hits.
- ✅ `cargo check -p cfd-math` passes.
- ✅ `cargo nextest run -p cfd-math interpolation` passes.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in `multigrid/interpolation.rs`.
- [x] Confirm the bounded slice preserves current nalgebra sparse/vector APIs.

### Phase 2: Execution (10-50%)
- [x] Replace interpolation scalar constants and index-distance conversion with Eunomia.
- [x] Replace interpolation quality metric extraction and absolute-value paths with Eunomia.
- [x] Remove touched test debug output and direct index-distance casts.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `multigrid/interpolation.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math interpolation` (14/14 passed).

---

# Sprint 1.96.71 Checklist: cfd-math Eunomia Multigrid Coarsening
**Goal**: Remove direct scalar conversion provider residue from multigrid
coarsening algorithms and quality analysis.

**Success Criteria**:
- ✅ Coarsening algorithms use Eunomia `FloatElement` for scalar constants and
  count-to-scalar conversion.
- ✅ Coarsening quality analysis uses Eunomia `FloatElement`/`NumericElement`
  for sentinels, distance thresholds, absolute values, and f64 metric
  extraction.
- ✅ Strength-matrix construction routes absolute-value dispatch through
  Eunomia.
- ✅ Value-semantic tests cover expected strength-matrix grid connectivity.
- ✅ `multigrid/coarsening` has no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, direct `T::from_f64`, direct `T::from_usize`, conversion
  fallback, `from_f64_or`, `SafeFromF64`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in `multigrid/coarsening`.
- [x] Confirm the bounded slice preserves current nalgebra sparse/vector APIs.

### Phase 2: Execution (10-50%)
- [x] Replace coarsening algorithm scalar constants and count conversions with Eunomia.
- [x] Replace coarsening quality sentinel, threshold, absolute-value, and f64 extraction paths with Eunomia.
- [x] Add value-semantic strength-matrix connectivity coverage.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `multigrid/coarsening`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math coarsening` (10/10 passed).

---

# Sprint 1.96.70 Checklist: cfd-math Eunomia Multigrid Smoothers
**Goal**: Remove direct scalar conversion fallbacks from multigrid smoother
owner paths.

**Success Criteria**:
- ✅ Multigrid smoother thresholds use Eunomia scalar conversion and
  absolute-value dispatch.
- ✅ Chebyshev eigenvalue defaults and recurrence constants use Eunomia scalar
  conversion.
- ✅ Immediate AMG smoother-owner threshold, relaxation, and complexity filters
  use Eunomia helpers.
- ✅ Value-semantic tests cover Gauss-Seidel, Jacobi, symmetric Gauss-Seidel,
  SOR, and Chebyshev smoother behavior.
- ✅ `smoothers.rs` has no direct `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`,
  `SafeFromF64`, `from_f64_or`, fallback, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in multigrid smoother owner paths.
- [x] Confirm the bounded slice preserves current nalgebra sparse/vector APIs.

### Phase 2: Execution (10-50%)
- [x] Replace smoother diagonal thresholds with Eunomia.
- [x] Replace Chebyshev eigenvalue/default constants with Eunomia.
- [x] Replace immediate AMG smoother-owner constants and complexity filters with Eunomia.
- [x] Add value-semantic smoother update/eigenvalue tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scans over `smoothers.rs` and `amg.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math test_gauss_seidel_smoother test_jacobi_smoother test_symmetric_gauss_seidel test_sor_smoother test_chebyshev_smoother` (5/5 passed).

---

# Sprint 1.96.69 Checklist: cfd-math Eunomia Convergence Monitor
**Goal**: Remove direct scalar conversion fallbacks from linear-solver convergence helpers.

**Success Criteria**:
- ✅ `ConvergenceMonitor::convergence_factor` uses Eunomia scalar conversion and `FloatElement::powf`.
- ✅ `ConvergenceMonitor::cg_theoretical_bound` uses Eunomia scalar conversion.
- ✅ `ConvergenceMonitor::validate_convergence` uses Eunomia scalar conversion.
- ✅ Value-semantic tests cover geometric convergence factor, CG theoretical bound, and validation rejection.
- ✅ `linear_solver/traits.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, `SafeFromF64`, `from_f64_or`, fallback, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in `linear_solver/traits.rs`.
- [x] Confirm the bounded slice preserves the current nalgebra vector/operator trait API.

### Phase 2: Execution (10-50%)
- [x] Replace convergence-factor exponent construction with Eunomia.
- [x] Replace CG theoretical-bound factor construction with Eunomia.
- [x] Replace validation safety-multiplier construction with Eunomia.
- [x] Add value-semantic convergence monitor tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `traits.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math convergence_factor cg_theoretical_bound validate_convergence` (4/4 passed, including the existing AMG convergence-factor match).

---

# Sprint 1.96.68 Checklist: cfd-math Eunomia Linear Operators
**Goal**: Remove direct `num-traits` scalar conversion from CPU finite-difference linear operators.

**Success Criteria**:
- ✅ Poisson/Laplacian operator constants use Eunomia `FloatElement`.
- ✅ Momentum/Energy operator constants use Eunomia `FloatElement`.
- ✅ Value-semantic center-stencil tests cover 2D Laplacian, 3D Poisson, 1D momentum, 2D momentum, and 2D energy operators.
- ✅ `operators/{poisson.rs,momentum.rs}` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in `linear_solver/operators/{poisson.rs,momentum.rs}`.
- [x] Confirm the bounded slice preserves the current nalgebra vector operator API.

### Phase 2: Execution (10-50%)
- [x] Replace Poisson/Laplacian scalar constants and provider bounds with Eunomia.
- [x] Replace momentum/energy scalar constants and provider bounds with Eunomia.
- [x] Add value-semantic finite-difference operator tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `poisson.rs` and `momentum.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math laplacian_center_impulse poisson_center_impulse` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-math momentum_1d_applies momentum_2d_applies energy_2d_applies` (3/3 passed).

---

# Sprint 1.96.67 Checklist: cfd-math Eunomia Schwarz/Cholesky Preconditioners
**Goal**: Remove direct `num-traits` scalar/provider residue from Schwarz and IncompleteCholesky preconditioners.

**Success Criteria**:
- ✅ Schwarz no longer imports or requires direct `num_traits::FromPrimitive`.
- ✅ Cholesky symmetry tolerance construction uses Eunomia `FloatElement`.
- ✅ Cholesky symmetry residual absolute-value dispatch uses Eunomia `NumericElement`.
- ✅ `schwarz.rs` and `cholesky.rs` have no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct provider residue in `linear_solver/preconditioners/{schwarz.rs,cholesky.rs}`.
- [x] Confirm the bounded slice preserves the current nalgebra CSR/vector preconditioner API.

### Phase 2: Execution (10-50%)
- [x] Remove Schwarz's stale direct `FromPrimitive` import and bound.
- [x] Replace Cholesky symmetry tolerance construction with Eunomia.
- [x] Route Cholesky symmetry residual absolute-value dispatch through Eunomia.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `schwarz.rs` and `cholesky.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math cholesky` (4/4 passed).
- [x] Run `cargo nextest run -p cfd-math schwarz` (1/1 passed).

---

# Sprint 1.96.66 Checklist: cfd-math Eunomia SSOR Preconditioner
**Goal**: Remove direct `num-traits` scalar conversion from the SSOR preconditioner.

**Success Criteria**:
- ✅ SSOR default relaxation construction uses Eunomia `FloatElement`.
- ✅ SSOR preconditioner impl bounds no longer require direct `num_traits::FromPrimitive`.
- ✅ `ssor.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct conversion residue in `linear_solver/preconditioners/ssor.rs`.
- [x] Confirm the bounded slice preserves the current nalgebra CSR/vector preconditioner API.

### Phase 2: Execution (10-50%)
- [x] Replace SSOR default relaxation construction.
- [x] Replace SSOR provider bounds with Eunomia `FloatElement`.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `ssor.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math ssor` (4/4 passed).

---

# Sprint 1.96.65 Checklist: cfd-math Eunomia Basic Preconditioners
**Goal**: Remove direct `num-traits` scalar conversions from basic Jacobi/SOR preconditioners.

**Success Criteria**:
- ✅ Jacobi diagonal tolerance construction uses Eunomia `FloatElement`.
- ✅ SOR omega bound and 1D-Poisson omega construction use Eunomia `FloatElement`.
- ✅ Basic preconditioner absolute-value dispatch uses Eunomia `NumericElement`.
- ✅ `basic.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, fallback, `from_f64_or`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct conversion residue in `linear_solver/preconditioners/basic.rs`.
- [x] Confirm the bounded slice preserves the current nalgebra CSR/vector preconditioner API.

### Phase 2: Execution (10-50%)
- [x] Replace Jacobi diagonal tolerance construction.
- [x] Replace SOR omega constants and omega construction.
- [x] Replace preconditioner absolute-value dispatch with Eunomia.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over `basic.rs`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math jacobi` (5/5 passed).
- [x] Run `cargo nextest run -p cfd-math sor` (6/6 passed).

---

# Sprint 1.96.64 Checklist: cfd-math Eunomia GMRES Chain
**Goal**: Remove direct `num-traits` scalar bounds from the GMRES-centered linear-solver chain.

**Success Criteria**:
- ✅ GMRES solver and Arnoldi helper use Eunomia `FloatElement`/`NumericElement`.
- ✅ `LinearSolverChain`, `DirectSparseSolver`, and block/SIMPLE preconditioners use Eunomia provider bounds.
- ✅ The touched slice has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `Float::`, direct `T::from_f64`, `T::from_usize`, conversion fallback, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit the GMRES-centered linear-solver call path.
- [x] Confirm the bounded slice preserves current nalgebra and rsparse backend surfaces.

### Phase 2: Execution (10-50%)
- [x] Replace GMRES and Arnoldi scalar bounds.
- [x] Replace linear-solver chain scalar bounds.
- [x] Replace direct sparse solver conversion bounds with Eunomia conversion dispatch.
- [x] Replace block/SIMPLE preconditioner scalar safeguards with Eunomia helpers.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run focused residue scan over the touched slice.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math gmres` (21/21 passed).
- [x] Run `cargo nextest run -p cfd-math direct_solver` (3/3 passed).
- [x] Run `cargo nextest run -p cfd-math block_preconditioner` (2/2 passed).

---

# Sprint 1.96.63 Checklist: cfd-math Eunomia Linear-Solver Config
**Goal**: Remove direct `num-traits` scalar conversion from iterative linear-solver default configuration.

**Success Criteria**:
- ✅ `IterativeSolverConfig::default` uses Eunomia `FloatElement` for tolerance construction.
- ✅ CG, BiCGSTAB, and GMRES default-construction paths declare the Atlas scalar provider bound.
- ✅ Linear-solver config docs name Moirai instead of Rayon for parallel SpMV.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit the bounded default-config residue in `linear_solver`.
- [x] Confirm the slice preserves the current nalgebra linear-solver API.

### Phase 2: Execution (10-50%)
- [x] Replace config default tolerance construction.
- [x] Propagate default-construction provider bounds through CG, BiCGSTAB, and GMRES.
- [x] Correct the stale parallel SpMV provider wording.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math default_solver` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-math test_gmres_configurable_trait` (1/1 passed).
- [x] Run focused source scan for direct conversion and stale Rayon residue.

---

# Sprint 1.96.62 Checklist: cfd-math Eunomia SIMD Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math SIMD utilities.

**Success Criteria**:
- ✅ CFD SIMD central-difference constants use Eunomia `FloatElement`.
- ✅ CFD SIMD field norm square-root dispatch uses Eunomia `NumericElement`.
- ✅ SIMD source has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Signed`, `Float::`, direct `T::from_f64`, conversion-fallback,
  `T::epsilon`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `simd`.
- [x] Confirm the bounded slice preserves the current Moirai-backed execution
  and nalgebra scalar API.

### Phase 2: Execution (10-50%)
- [x] Replace CFD SIMD central-difference constants.
- [x] Replace CFD field norm square-root dispatch.
- [x] Remove direct `num_traits::FromPrimitive` bounds from SIMD CFD helpers.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math simd` (26/26 passed).
- [x] Run SIMD source scan for direct conversion and stale Rayon residue.

---

# Sprint 1.96.61 Checklist: cfd-math Eunomia Sparse Scalars
**Goal**: Remove direct `num-traits` scalar conversions and math dispatch from cfd-math sparse utilities.

**Success Criteria**:
- ✅ Sparse stencil constants use Eunomia `FloatElement`.
- ✅ Frobenius norm, condition-estimate absolute values, diagonal dominance,
  and singular-diagonal thresholds dispatch through Eunomia.
- ✅ Sparse source has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `Signed`, `Float::`, direct `T::from_f64`, conversion-fallback,
  `T::epsilon`, or stale `rayon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `sparse`.
- [x] Confirm the bounded slice preserves the current CSR/DVector sparse API.

### Phase 2: Execution (10-50%)
- [x] Replace sparse stencil constants.
- [x] Replace Frobenius norm and condition-estimate scalar math.
- [x] Replace diagonal dominance absolute-value dispatch.
- [x] Correct sparse SpMV docs to name Moirai instead of Rayon.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math sparse` (15/15 passed).
- [x] Run sparse source scan for direct conversion and stale Rayon residue.

---

# Sprint 1.96.60 Checklist: cfd-math Eunomia Nonlinear Solver Scalars
**Goal**: Remove direct `num-traits` scalar conversions and math dispatch from cfd-math Anderson/JFNK nonlinear solvers.

**Success Criteria**:
- ✅ Anderson/JFNK default constants use Eunomia `FloatElement`.
- ✅ QR/JFNK square-root, absolute-value, min/max clamp, and finite-difference
  perturbation math dispatch through Eunomia `NumericElement`/`FloatElement`.
- ✅ Nonlinear-solver source has no direct `num_traits`, `FromPrimitive`,
  `Float::`, direct `T::from_f64`, conversion-fallback, or `T::epsilon` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `nonlinear_solver`.
- [x] Record that the scalar-only slice preserved nalgebra vector ownership at
  the time; superseded by Sprint 1.96.109 Leto vector migration.

### Phase 2: Execution (10-50%)
- [x] Add Eunomia scalar helpers for Anderson/JFNK defaults.
- [x] Replace Anderson QR norm and diagonal checks.
- [x] Replace JFNK perturbation, EW forcing, Givens, and back-substitution math.

### Phase 3: Closure (50%+)
- [x] Run touched-file rustfmt.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math nonlinear_solver` (9/9 passed).
- [x] Run nonlinear-solver source scan for direct conversion residue.

---

# Sprint 1.96.59 Checklist: cfd-math Eunomia Pressure-Velocity Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math SIMPLE pressure-velocity defaults.

**Success Criteria**:
- ✅ SIMPLE default tolerance and relaxation constants use Eunomia
  `FloatElement`.
- ✅ `SIMPLEConfig::new` no longer requires a scalar-conversion trait.
- ✅ Pressure-velocity source has no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_f64`, `T::from_usize`, `T::from_f32`,
  conversion-fallback, or conversion `expect` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `pressure_velocity`.
- [x] Confirm the bounded slice preserves the current SIMPLE config/result API.

### Phase 2: Execution (10-50%)
- [x] Add a Eunomia scalar helper for SIMPLE defaults.
- [x] Replace SIMPLE default constants.
- [x] Narrow explicit SIMPLE config construction bounds.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math pressure_velocity` (2/2 passed).
- [x] Run pressure-velocity source scan for direct conversion residue.

---

# Sprint 1.96.58 Checklist: cfd-math Eunomia Iterator Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math iterator utilities.

**Success Criteria**:
- ✅ Stencil coefficients use Eunomia `FloatElement`.
- ✅ Iterator statistics count conversion uses Eunomia `FloatElement` and
  `NumericElement`.
- ✅ Standard-deviation square-root dispatch is explicitly routed through
  Eunomia.
- ✅ Every declared 3-point stencil pattern returns real second-derivative
  coefficients instead of a placeholder zero vector.
- ✅ Iterator source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `iterators`.
- [x] Confirm the bounded slice preserves the current iterator trait API.

### Phase 2: Execution (10-50%)
- [x] Add Eunomia scalar helpers for stencil constants and iterator counts.
- [x] Replace stencil coefficient scalar construction.
- [x] Replace iterator statistics scalar conversion and math dispatch.
- [x] Replace the unsupported second-derivative zero branch with real 3-point
  coefficients.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math iterators` (7/7 passed).
- [x] Run iterator source scan for direct conversion residue.

---

# Sprint 1.96.57 Checklist: cfd-math Eunomia WENO Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math WENO high-order reconstruction.

**Success Criteria**:
- ✅ WENO5/WENO7 epsilon defaults use Eunomia `FloatElement`.
- ✅ WENO5/WENO7 linear weights, ENO coefficients, and smoothness constants use
  Eunomia `FloatElement`.
- ✅ WENO denominator squaring avoids the `FloatElement`/nalgebra `powi`
  ambiguity.
- ✅ High-order source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `high_order`.
- [x] Confirm the bounded slice preserves the current WENO reconstruction API.

### Phase 2: Execution (10-50%)
- [x] Add shared WENO Eunomia scalar and square helpers.
- [x] Replace WENO5 scalar constants and `FromPrimitive` bounds.
- [x] Replace WENO7 scalar constants and `FromPrimitive` bounds.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math weno` (6/6 passed).
- [x] Run high-order source scan for direct conversion residue.

---

# Sprint 1.96.56 Checklist: cfd-math Eunomia Interpolation Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math cubic-spline interpolation.

**Success Criteria**:
- ✅ Cubic-spline Thomas-algorithm constants use Eunomia `FloatElement`.
- ✅ Cubic-spline trait bounds no longer require `num_traits::FromPrimitive`.
- ✅ Cubic-spline source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `T::from_f32`, `From<f64>`, or
  conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `interpolation`.
- [x] Confirm the bounded slice preserves the current interpolation trait API.

### Phase 2: Execution (10-50%)
- [x] Add a Eunomia scalar helper for cubic-spline constants.
- [x] Replace Thomas-algorithm scalar constants.
- [x] Replace cubic-spline `FromPrimitive` bounds with `FloatElement`.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math interpolation` (14/14 passed).
- [x] Run cubic-spline source scan for direct conversion residue.

---

# Sprint 1.96.55 Checklist: cfd-math Eunomia Differentiation Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math differentiation operators.

**Success Criteria**:
- ✅ Finite-difference stencil constants use Eunomia `FloatElement`.
- ✅ Gradient/divergence/curl constants use Eunomia `FloatElement`.
- ✅ SIMD helper scalar staging uses Eunomia rather than `T::from_f32`.
- ✅ Differentiation source has no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_f64`, `T::from_usize`, `From<f64>`, `T::from_f32`,
  or conversion-fallback hits in the touched operator files.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `differentiation`.
- [x] Confirm the bounded slice preserves the current nalgebra
  `DVector`/`Vector3` differentiation API.

### Phase 2: Execution (10-50%)
- [x] Add Eunomia scalar helpers for differentiation constants.
- [x] Replace finite-difference stencil constants.
- [x] Replace gradient/divergence/curl constants.
- [x] Replace SIMD helper scalar staging with Eunomia conversions.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math differentiation` (12/12 passed).
- [x] Run touched differentiation source scan for direct conversion residue.

---

# Sprint 1.96.54 Checklist: cfd-math Eunomia Integration Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math integration quadrature.

**Success Criteria**:
- ✅ Trapezoidal, Simpson, and Gauss quadrature constants use Eunomia
  `FloatElement`.
- ✅ Composite interval index/count conversions use Eunomia helpers.
- ✅ Adaptive tolerance/error dispatch uses Eunomia
  `FloatElement`/`NumericElement`.
- ✅ Tetrahedral quadrature constants use Eunomia.
- ✅ Integration source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, `From<f64>`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `integration`.
- [x] Confirm the bounded slice preserves the current `Quadrature` trait and
  nalgebra `RealField` API.

### Phase 2: Execution (10-50%)
- [x] Add Eunomia scalar helpers for integration constants and exact counts.
- [x] Replace 1D quadrature constants and weights.
- [x] Replace composite/adaptive/tensor scalar conversion bounds.
- [x] Replace tetrahedral quadrature constants.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math integration` (5/5 passed).
- [x] Run integration source scan for direct conversion residue.

---

# Sprint 1.96.53 Checklist: cfd-math Eunomia Exponential Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math exponential time stepping.

**Success Criteria**:
- ✅ ERK4 stage weights and combination denominator use Eunomia `FloatElement`.
- ✅ phi-function small-argument thresholds use Eunomia `FloatElement`.
- ✅ scaling/squaring factorial conversion uses a Eunomia helper.
- ✅ Exponential source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::exponential`.
- [x] Confirm the bounded slice preserves the current `DVector`/`DMatrix` state API.

### Phase 2: Execution (10-50%)
- [x] Add local Eunomia scalar helpers.
- [x] Replace ERK4 scalar coefficients.
- [x] Replace phi-function small-argument thresholds.
- [x] Replace scaling/squaring factorial conversion.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math time_stepping::exponential::` (2/2 passed).
- [x] Run touched source scan for direct conversion residue.

---

# Sprint 1.96.52 Checklist: cfd-math Eunomia IMEX Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math IMEX time stepping.

**Success Criteria**:
- ✅ Newton tolerance uses Eunomia `FloatElement`.
- ✅ ARS343 gamma/delta and tableau coefficients use Eunomia helpers.
- ✅ IMEX source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::imex`.
- [x] Confirm the bounded slice preserves the current `DVector`/`DMatrix` state API.

### Phase 2: Execution (10-50%)
- [x] Add local Eunomia scalar helper.
- [x] Replace Newton tolerance conversion.
- [x] Replace ARS343 gamma/delta constants and scalar square-root dispatch.
- [x] Replace explicit/implicit tableau coefficients and solution weights.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math imex` (5/5 passed).
- [x] Run touched source scan for direct conversion residue.

---

# Sprint 1.96.51 Checklist: cfd-math Eunomia RKC Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math Runge-Kutta-Chebyshev time stepping.

**Success Criteria**:
- ✅ RKC defaults use Eunomia `FloatElement`.
- ✅ Chebyshev recurrence constants and count conversions use Eunomia helpers.
- ✅ Adaptive error-control constants and power dispatch use Eunomia.
- ✅ RKC source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, `T::from_usize`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::rk_chebyshev`.
- [x] Confirm the bounded slice preserves the current `DVector` state API.

### Phase 2: Execution (10-50%)
- [x] Add local Eunomia scalar helpers for floating constants and exact count conversion.
- [x] Replace RKC default tolerance, damping, and safety constants.
- [x] Replace Chebyshev recurrence constants and stage/vector length conversions.
- [x] Route `sqrt`, `ln`, `sinh`, `cosh`, `abs`, and `powf` dispatch through Eunomia.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math rk_chebyshev` (4/4 passed).
- [x] Run touched source scan for direct conversion residue.

---

# Sprint 1.96.50 Checklist: cfd-math Eunomia Adaptive Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math adaptive time stepping.

**Success Criteria**:
- ✅ Adaptive defaults use Eunomia `FloatElement`.
- ✅ PI controller constants and clamp bounds use Eunomia `FloatElement`.
- ✅ Dormand-Prince tableau coefficients use Eunomia `FloatElement`.
- ✅ Adaptive source has no `num_traits`, `FromPrimitive`, `ToPrimitive`,
  `T::from_f64`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::adaptive`.
- [x] Confirm the bounded slice preserves the current `DVector` state API.

### Phase 2: Execution (10-50%)
- [x] Add a local Eunomia `FloatElement` scalar helper.
- [x] Replace adaptive-stepper defaults and rejection scaling.
- [x] Replace PI controller gains, power dispatch, and clamp bounds.
- [x] Replace Dormand-Prince c/a/b4/b5 tableau constants.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math adaptive` (3/3 passed).
- [x] Run touched source scan for direct conversion residue.

---

# Sprint 1.96.49 Checklist: cfd-math Eunomia Runge-Kutta Scalars
**Goal**: Remove direct `num-traits` scalar conversions from cfd-math Runge-Kutta methods and repair low-storage RK4 value semantics.

**Success Criteria**:
- ✅ RK3/RK4/low-storage RK4 constants use Eunomia `FloatElement`.
- ✅ `LowStorageRK4` preserves constant solutions for zero RHS.
- ✅ Runge-Kutta source has no `num_traits`, `FromPrimitive`,
  `ToPrimitive`, `T::from_f64`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::runge_kutta`.
- [x] Identify low-storage RK4 recurrence as a touched correctness defect.

### Phase 2: Execution (10-50%)
- [x] Add a local Eunomia `FloatElement` scalar helper.
- [x] Replace RK3, RK4, and Carpenter-Kennedy coefficient conversions.
- [x] Implement the Carpenter-Kennedy 2N residual recurrence.
- [x] Replace weak low-storage RK4 assertions with value-semantic tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math runge_kutta` (5/5 passed).
- [x] Run touched source scan for direct conversion residue.

---

# Sprint 1.96.48 Checklist: cfd-math Eunomia Stability Scalars
**Goal**: Remove direct `num-traits` scalar conversions from the cfd-math stability analyzer.

**Success Criteria**:
- ✅ Stability analyzer constants use Eunomia `FloatElement`.
- ✅ Stability diagnostic numeric formatting uses Eunomia `NumericElement`.
- ✅ Stability source files have no `num_traits`, `ToPrimitive`,
  `FromPrimitive`, `T::from_f64`, or conversion-fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit scalar conversion residue in `time_stepping::stability`.
- [x] Confirm the bounded slice does not attempt the larger nalgebra/Leto
  Butcher-tableau migration.

### Phase 2: Execution (10-50%)
- [x] Add local Eunomia conversion helpers over `FloatElement` and
  `NumericElement`.
- [x] Replace CFL, RK-order, von-Neumann, and recommendation scalar conversions.
- [x] Disambiguate numeric `abs` dispatch through Eunomia.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math stability` (5/5 passed).
- [x] Run touched stability source scan for direct conversion residue.
- [x] Run `git diff --check`.
- [x] Record package-level cfd-math fmt residual outside the touched files.

---

# Sprint 1.96.46 Checklist: cfd-2d Eunomia Scheme Complex
**Goal**: Remove direct `num-complex` ownership from the cfd-2d scheme boundary.

**Success Criteria**:
- ✅ Scheme amplification factors use `eunomia::Complex<f64>`.
- ✅ `crates/cfd-2d` has no direct `num_complex`, `num-complex`, or `NumComplex` source/manifest hits.
- ✅ cfd-2d provider-bound compile errors from cfd-core `FloatElement` requirements are closed.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct `num_complex` use in cfd-2d.
- [x] Confirm `SpatialDiscretization::amplification_factor` is the bounded complex-provider API.
- [x] Identify the cfd-core `FloatElement` requirements blocking cfd-2d compilation.

### Phase 2: Execution (10-50%)
- [x] Add direct Eunomia dependency to cfd-2d.
- [x] Replace the scheme amplification factor return and constructors with `eunomia::Complex`.
- [x] Remove direct `num-complex` manifest ownership from cfd-2d.
- [x] Add explicit `FloatElement` bounds where 2D paths construct cfd-core grids or solver configs.

### Phase 3: Closure (50%+)
- [x] Run `cargo fmt --package cfd-2d --package cfd-3d --package cfd-validation --check`.
- [x] Run `cargo check -p cfd-2d`.
- [x] Run `cargo check -p cfd-3d`.
- [x] Run `cargo check -p cfd-validation`.
- [x] Run `cargo nextest run -p cfd-2d` (563/563 passed, 27 skipped).
- [x] Run cfd-2d source/manifest scan for `num_complex`, `num-complex`, and `NumComplex`.
- [x] Record residual direct `num-traits`/nalgebra provider migration risk.

---

# Sprint 1.96.45 Checklist: cfd-math/cfd-validation Eunomia Stability Complex
**Goal**: Remove direct `num-complex` ownership from the time-stepping stability-analysis boundary.

**Success Criteria**:
- ✅ Von Neumann spatial-operator callbacks use `eunomia::Complex<f64>`.
- ✅ Explicit RK stability evaluation avoids `num_complex` dense matrix allocation/inversion.
- ✅ `crates/cfd-math` and `crates/cfd-validation` have no `num_complex`, `num-complex`, or `NumComplex` source/manifest hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct `num_complex` use in cfd-math and cfd-validation.
- [x] Confirm the existing analyzer accepts only explicit lower-triangular RK tableaux.
- [x] Select the stability analyzer and validation stability closures as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Add direct Eunomia dependencies to cfd-math and cfd-validation.
- [x] Replace `num_complex::Complex<f64>` callback types and constructors with `eunomia::Complex<f64>`.
- [x] Replace the dense complex matrix inverse with forward substitution over `(I - zA)x = 1`.
- [x] Remove direct `num-complex` manifest dependencies from the touched crates.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-math`.
- [x] Run `cargo nextest run -p cfd-math stability` (5/5 passed).
- [x] Run touched-crate source/manifest scan.
- [x] Record that the downstream cfd-validation blocker was cleared in the later cfd-2d/cfd-validation provider-bound slice.

---

# Sprint 1.96.44 Checklist: cfd-core Eunomia Rhie-Chow Interpolation
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core Rhie-Chow interpolation.

**Success Criteria**:
- ✅ Default relaxation and two-face constants convert through Eunomia `FloatElement`.
- ✅ `physics/fluid_dynamics/rhie_chow.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion fallback hits.
- ✅ U-face and V-face pressure-correction formulas have value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit Rhie-Chow scalar conversion sites.
- [x] Confirm the bounded slice preserves existing pressure-correction formulas.
- [x] Select `physics/fluid_dynamics/rhie_chow.rs` as the Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and bounds with Eunomia `FloatElement`.
- [x] Route relaxation and two constants through Eunomia.
- [x] Add value-semantic u-face and v-face interpolation tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features rhie_chow` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (176/176 passed).
- [x] Run touched-file source scan.

---

# Sprint 1.96.43 Checklist: cfd-core Eunomia Boundary Geometry
**Goal**: Remove direct scalar conversion fallbacks from cfd-core boundary geometry measure formulas.

**Success Criteria**:
- ✅ Sphere and cylinder measure constants convert through Eunomia `FloatElement`.
- ✅ `contains_point` and `dimension` keep their existing conversion-free scalar bounds.
- ✅ `physics/boundary/geometry.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion fallback hits.
- ✅ Line, sphere, cylinder, and unsupported-measure behavior has value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit boundary geometry scalar conversion sites.
- [x] Select `physics/boundary/geometry.rs` as the bounded Eunomia migration slice.
- [x] Preserve existing formulas and avoid widening conversion-free method bounds.

### Phase 2: Execution (10-50%)
- [x] Route sphere and cylinder constants through Eunomia `FloatElement`.
- [x] Split the `BoundaryGeometry` impl so only `measure()` carries the provider conversion bound.
- [x] Add value-semantic line, sphere, cylinder, and unsupported-measure tests.

### Phase 3: Closure (50%+)
- [x] Run touched-file `rustfmt --check`.
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features boundary::geometry` (4/4 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (174/174 passed).
- [x] Run touched-file source scan.

---

# Sprint 1.96.42 Checklist: cfd-core Eunomia Boundary Ghost Cells
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core boundary ghost-cell formulas.

**Success Criteria**:
- ✅ Ghost-cell Dirichlet, Neumann, and Robin constants convert through Eunomia `FloatElement`.
- ✅ Robin singularity reporting uses Eunomia `NumericElement` instead of `num_traits::ToPrimitive`.
- ✅ Degenerate Robin coefficients reject with `BoundaryErrorKind::RobinSingularity` before division.
- ✅ `physics/boundary/ghost_cells.rs` has no direct `num_traits`, `FromPrimitive`, `ToPrimitive`, `T::from_`, or conversion fallback hits.
- ✅ Second-order, fourth-order, and degenerate Robin ghost-cell behavior has value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit boundary ghost-cell scalar conversion and error-reporting sites.
- [x] Confirm the bounded slice preserves existing ghost-cell formulas.
- [x] Select `physics/boundary/ghost_cells.rs` as the Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive`/`ToPrimitive` imports and bounds with Eunomia `FloatElement`/`NumericElement`.
- [x] Route Dirichlet, Neumann, Robin, and singularity tolerance constants through Eunomia.
- [x] Return a typed Robin singularity error when both Robin coefficients are zero.
- [x] Add value-semantic fourth-order Dirichlet, fourth-order Neumann, and degenerate Robin tests.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features ghost_cells` (5/5 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (170/170 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.41 Checklist: cfd-core Eunomia Staggered Grid
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core staggered-grid coordinate construction.

**Success Criteria**:
- ✅ Staggered-grid dimensions and coordinate indices convert through Eunomia `FloatElement`.
- ✅ Integer grid dimensions and indices assert exact representability before scalar conversion.
- ✅ `geometry/staggered.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_`, `from_usize`, or conversion fallback hits.
- ✅ Uniform and stretched staggered-grid coordinate behavior remains value-semantic.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit staggered-grid scalar conversion sites.
- [x] Confirm the bounded slice preserves existing uniform and stretched grid formulas.
- [x] Select `geometry/staggered.rs` as the Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and bounds with Eunomia `FloatElement`.
- [x] Add an exact-representability helper for integer grid dimensions and indices.
- [x] Route half constants through Eunomia `FloatElement`.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features staggered` (5/5 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (167/167 passed).
- [x] Run touched-file `rustfmt --check`, `git diff --check`, and source scan.

---

# Sprint 1.96.40 Checklist: cfd-core Eunomia Boundary Time Functions
**Goal**: Remove direct `num-traits` scalar conversion and nalgebra math dispatch from cfd-core boundary time functions and ghost-cell constants.

**Success Criteria**:
- ✅ Boundary time functions use Eunomia `FloatElement` for `2π`, `sin`, and `exp`.
- ✅ Boundary ghost-cell reflection uses Eunomia `FloatElement` for the `2` constant.
- ✅ Boundary applicator, concrete applicator, specification, and manager surfaces carry the provider bound required by time-dependent evaluation.
- ✅ Touched boundary files have no direct `num_traits`, `FromPrimitive`, `T::from_f64`, conversion fallback, bare `.exp()`, or bare `.sin()` hits.
- ✅ Sinusoidal, exponential, boundary scaling, and ghost-cell behavior have value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit boundary time-function and ghost-cell conversion sites.
- [x] Confirm the bounded slice preserves existing formulas and public boundary behavior.
- [x] Select `physics/boundary/{time_dependent,applicator,applicators,manager,specification}.rs` as the Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and fallback conversions with Eunomia `FloatElement`.
- [x] Route sinusoidal and exponential time-function math through Eunomia.
- [x] Propagate `FloatElement` bounds through boundary applicator/specification/manager surfaces.
- [x] Add value-semantic tests for sinusoidal, exponential, boundary scaling, and ghost-cell behavior.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features boundary` (30/30 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (167/167 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.39 Checklist: cfd-core Eunomia Temperature-Dependent Fluids
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core temperature-dependent fluid models.

**Success Criteria**:
- ✅ Polynomial fluid models no longer carry stale scalar-conversion bounds.
- ✅ Arrhenius, Andrade, and Sutherland math dispatch uses Eunomia `FloatElement`.
- ✅ `physics/fluid/temperature.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_f64`, conversion fallback, bare `.exp()`, or bare `.powf(...)` hits.
- ✅ Polynomial, Andrade, Arrhenius, and Sutherland behavior has value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit temperature-dependent fluid conversion bounds and scalar literal sites.
- [x] Confirm the bounded slice preserves existing fluid APIs and formulas.
- [x] Select `physics/fluid/temperature.rs` as the next Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Remove stale `FromPrimitive` import and bounds.
- [x] Route Arrhenius/Andrade exponentials through Eunomia `FloatElement`.
- [x] Route Sutherland exponent conversion and `powf` through Eunomia `FloatElement`.
- [x] Add value-semantic tests for polynomial density/viscosity, Andrade viscosity/domain, and Sutherland formula evaluation.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features fluid::temperature` (6/6 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (162/162 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.38 Checklist: cfd-core Eunomia Hemolysis Constants
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core hemolysis calculator and platelet activation constants.

**Success Criteria**:
- ✅ NIH, MIH, exposure-time, and platelet activation constants use Eunomia `FloatElement`.
- ✅ `physics/hemolysis/{calculator,trauma}.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits.
- ✅ NIH, MIH, exposure time, and platelet activation probability have value-semantic tests.
- ✅ Platelet activation above threshold returns a probability in `(0, 1)` following `1 - exp(-k · excess_stress · exposure_time)`.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit hemolysis calculator and platelet activation conversion sites.
- [x] Confirm the bounded slice changes scalar constants and the activation probability sign defect.
- [x] Select `physics/hemolysis/{calculator,trauma}.rs` as the next Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and bounds with Eunomia `FloatElement`.
- [x] Replace NIH, MIH, exposure-time, and platelet activation constants with Eunomia conversions.
- [x] Correct platelet activation probability to use the decaying exponential.
- [x] Add value-semantic tests for NIH, MIH, exposure time, and platelet activation probability.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features hemolysis` (9/9 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (160/160 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.37 Checklist: cfd-core Eunomia Mesh Quality Thresholds
**Goal**: Remove direct `num-traits` scalar conversion from cfd-core mesh quality thresholds.

**Success Criteria**:
- ✅ Aspect-ratio, skewness, and orthogonality threshold constants use Eunomia `FloatElement`.
- ✅ `geometry/mesh/quality.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits.
- ✅ Mesh quality level boundaries and recommendations have value-semantic tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit mesh quality threshold conversion sites and existing geometry coverage.
- [x] Confirm the bounded slice only changes fixed threshold constants and preserves comparison semantics.
- [x] Select `geometry/mesh/quality.rs` as the next Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and bounds with Eunomia `FloatElement`.
- [x] Replace aspect-ratio, skewness, and orthogonality threshold conversions with Eunomia constants.
- [x] Add value-semantic tests for best-quality, worst-quality, and strict-threshold boundary behavior.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features mesh_quality` (3/3 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (158/158 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.36 Checklist: cfd-core Eunomia CPU Backend
**Goal**: Remove direct `num-traits` scalar conversion from the cfd-core CPU backend.

**Success Criteria**:
- ✅ CPU advection domain-parameter conversion uses Eunomia `FloatElement`.
- ✅ The local `safe_f64_to_t` fallback helper is removed.
- ✅ `CpuBuffer` impls do not require scalar-conversion bounds.
- ✅ `compute/cpu.rs` has no direct `num_traits`, `FromPrimitive`, `safe_f64_to_t`, or conversion fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit CPU backend storage and advection conversion sites.
- [x] Confirm conversion is only required at the `KernelParams` f64 domain-parameter boundary.
- [x] Select `compute/cpu.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct `FromPrimitive` import and bounds with Eunomia `FloatElement` on the advection kernel.
- [x] Delete the fallback conversion helper.
- [x] Narrow `CpuBuffer` construction, buffer trait, and debug impl bounds to conversion-free scalar bounds.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features test_cpu_advection_kernel_linear_exactness` (1/1 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (155/155 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.35 Checklist: cfd-core Eunomia Time Integrators
**Goal**: Remove direct `num-traits` scalar conversions from cfd-core time integrators.

**Success Criteria**:
- ✅ RK2, RK4, and Crank-Nicolson runtime constants use Eunomia `FloatElement`.
- ✅ Backward Euler and Crank-Nicolson default tolerances use Eunomia constants with no silent `T::one()` fallback.
- ✅ `compute/time/integrators.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits.
- ✅ Explicit and implicit integrator behavior has value-semantic regression tests.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit time integrator scalar constants, default tolerances, and call sites.
- [x] Confirm `ForwardEuler` needs no scalar conversion seam.
- [x] Select `compute/time/integrators.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `FromPrimitive` bounds with Eunomia `FloatElement` where scalar constants are required.
- [x] Replace RK2/RK4/Crank-Nicolson constants with Eunomia conversions.
- [x] Replace implicit default tolerances with Eunomia constants.
- [x] Add value-semantic tests for explicit updates, defaults, invalid iteration configuration, and nonconvergence errors.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features compute::time::integrators` (6/6 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (155/155 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.34 Checklist: cfd-core Eunomia Time-Step Controllers
**Goal**: Remove direct `num-traits` scalar conversions and float helpers from cfd-core time-step controllers.

**Success Criteria**:
- ✅ Adaptive and variable controller defaults use Eunomia `FloatElement` constants.
- ✅ Controller runtime math uses Eunomia `FloatElement::powf` and local comparison helpers instead of `num_traits::Float`.
- ✅ Invalid integration order returns a typed configuration error instead of a silent fallback.
- ✅ `compute/time/controllers.rs` has no direct `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit time-step controller defaults, runtime math, and call sites.
- [x] Confirm there are no in-repo callers requiring call-site migration.
- [x] Select `compute/time/controllers.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `FromPrimitive`/`Float` bounds with Eunomia `FloatElement`.
- [x] Replace default constants and power functions with Eunomia calls.
- [x] Replace min/max helper calls with local value-preserving comparison helpers.
- [x] Return typed errors for invalid integration-order conversion.
- [x] Add value-semantic tests for defaults, invalid order rejection, and clamp behavior.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features compute::time::controllers` (4/4 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (149/149 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.33 Checklist: cfd-core Eunomia Solver Config Defaults
**Goal**: Remove direct `num-traits` scalar conversions from cfd-core solver configuration defaults.

**Success Criteria**:
- ✅ `SolverConfig::default` uses Eunomia `FloatElement` constants.
- ✅ `LinearSolverConfig::default` uses Eunomia `FloatElement` constants.
- ✅ `SolverConfigBuilder::new` delegates to the canonical `SolverConfig::default` implementation.
- ✅ `compute/solver/config.rs` has value-semantic tests for default/builder parity and no direct `num_traits` conversion hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit solver config default constants and builder duplication.
- [x] Confirm Eunomia `FloatElement::from_f64` covers all changed constants.
- [x] Select `compute/solver/config.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `FromPrimitive` with Eunomia `FloatElement` in solver config defaults/builders.
- [x] Replace default scalar conversions with Eunomia constants.
- [x] Consolidate `SolverConfigBuilder::new` through `SolverConfig::default`.
- [x] Add value-semantic tests for solver and linear-solver config defaults.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features solver_config` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (146/146 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.32 Checklist: cfd-core Eunomia Abstraction Defaults
**Goal**: Remove direct `num-traits` scalar conversions from cfd-core abstraction defaults.

**Success Criteria**:
- ✅ `FieldState::new` uses Eunomia `FloatElement` for its default time step.
- ✅ `ProblemParameters::default` uses Eunomia `FloatElement` for reference pressure and temperature.
- ✅ `abstractions/{state,problem}.rs` have no `num_traits`, `FromPrimitive`, `T::from_f64`, or zero/one conversion fallback hits.
- ✅ Existing abstraction tests and cfd-core no-default test suite remain green.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit abstraction default constants and builder construction.
- [x] Confirm Eunomia `FloatElement::from_f64` covers the default constants.
- [x] Select `abstractions/{state,problem}.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `FromPrimitive` with Eunomia `FloatElement` in constructor/default bounds.
- [x] Replace default scalar conversions with Eunomia constants.
- [x] Split `FieldState` impls so non-constructor field methods do not require scalar-conversion bounds.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features abstractions` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (144/144 passed).
- [x] Run touched-file `rustfmt --check` and source scan.

---

# Sprint 1.96.31 Checklist: cfd-core Eunomia Fluid Validation Thresholds
**Goal**: Remove direct `num-traits` threshold conversions from cfd-core fluid validation.

**Success Criteria**:
- ✅ `PropertyBounds::default` uses Eunomia `FloatElement` constants.
- ✅ Reynolds, Prandtl, temperature, and pressure validators use Eunomia threshold constants.
- ✅ `validation.rs` has no `num_traits`, `FromPrimitive`, `T::from_f64`, or conversion fallback hits.
- ✅ Existing validation tests and cfd-core no-default test suite remain green.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit fluid validation thresholds and tests.
- [x] Confirm Eunomia `FloatElement::from_f64` covers all validation constants.
- [x] Select `physics/fluid/validation.rs` as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `FromPrimitive` with Eunomia `FloatElement`.
- [x] Replace property-bound default conversions with Eunomia constants.
- [x] Replace Reynolds/Prandtl/temperature/pressure limits with Eunomia constants.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features physics::fluid::validation` (2/2 passed).
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (144/144 passed).
- [x] Run touched-file `rustfmt --check`, `git diff --check`, and source scan.

---

# Sprint 1.96.30 Checklist: cfd-core Eunomia Material/Fluid Constants
**Goal**: Remove direct `num-traits` scalar conversion patterns from cfd-core material constructors and their immediate fluid constructor dependencies.

**Success Criteria**:
- ✅ Material constructors use Eunomia `FloatElement` for scalar constants.
- ✅ Constant-property fluid/database constructors use Eunomia `FloatElement`.
- ✅ The touched cone has no `num_traits`, `FromPrimitive`, `T::from_f64`, or silent zero/one fallback conversion hits.
- ✅ Existing cfd-core no-default tests remain green.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit cfd-core material constructors, material database setup, and immediate fluid constructor dependencies.
- [x] Confirm Eunomia `FloatElement` and `NumericElement` cover the required constants and ideal-gas math helpers.
- [x] Select material/fluid constants as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace direct material `FromPrimitive` bounds with Eunomia `FloatElement`.
- [x] Replace fluid database conversion-error branches with Eunomia constants.
- [x] Replace constant-property fluid and ideal-gas conversion fallbacks with Eunomia constants/math.
- [x] Keep the existing database `Result` API while removing conversion-only failure branches.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo check -p cfd-core`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (144/144 passed).
- [x] Run touched-file `rustfmt --check` and `git diff --check`.
- [x] Run touched-cone source scan for removed direct `num-traits` conversion patterns.

---

# Sprint 1.96.29 Checklist: cfd-core Eunomia Physics Value Boundary
**Goal**: Remove direct `num-traits` scalar constant conversion patterns from cfd-core physics value objects and their immediate management aggregate callers.

**Success Criteria**:
- ✅ `cfd-core` physics value objects use Eunomia `FloatElement` for scalar constants.
- ✅ Touched management aggregates no longer use `FromPrimitive` or silent constant conversion fallbacks.
- ✅ Invalid custom Reynolds parameter construction is surfaced as `Result<Self>` instead of substituted with a default.
- ✅ The touched source cone has no `num_traits`, `FromPrimitive`, `T::from_f64`, or `unwrap_or_else(T::zero/one)` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit cfd-core physics value objects and dependent management aggregate call sites.
- [x] Confirm Eunomia exposes `FloatElement::from_f64` for scalar constants.
- [x] Select value objects plus immediate management aggregate callers as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Add direct cfd-core `eunomia.workspace = true` dependency.
- [x] Replace direct `FromPrimitive` bounds with Eunomia `FloatElement` in the touched cone.
- [x] Replace scalar constant conversion fallbacks with Eunomia conversions and explicit invariant expectations where construction is known-valid.
- [x] Make `PhysicalParameters::with_reynolds` return `Result<Self>` for invalid input instead of substituting a default.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --no-default-features`.
- [x] Run `cargo nextest run -p cfd-core --no-default-features` (144/144 passed).
- [x] Run touched-file `rustfmt --check`.
- [x] Run touched-cone source scan for removed `num_traits`/`FromPrimitive`/fallback conversion patterns.
- [x] Record residual cfd-core `num-traits`, nalgebra, and package-fmt drift for later slices.

---

# Sprint 1.96.28 Checklist: cfd-python Eunomia Blood Binding Boundary
**Goal**: Remove cfd-python's direct num-traits dependency by routing blood-model PyO3 conversions through Eunomia.

**Success Criteria**:
- ✅ `cfd-python` depends on `eunomia.workspace = true`, not `num-traits = "0.2"`.
- ✅ Blood-model PyO3 wrappers have no `FromPrimitive`/`ToPrimitive` use.
- ✅ Getter methods expose actual Rust model field values instead of hardcoded fallback defaults.
- ✅ `crates/cfd-python` has no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct cfd-python `num-traits` use sites in `blood.rs` and the manifest.
- [x] Confirm Eunomia exposes `NumericElement::to_f64` for primitive f64 fields.
- [x] Select blood-model bindings as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Add direct cfd-python `eunomia` dependency and remove direct `num-traits`.
- [x] Replace `FromPrimitive` shear-rate conversions with direct f64 model calls.
- [x] Replace getter `ToPrimitive` fallback conversions with Eunomia `NumericElement`.
- [x] Clear dependency-chain clippy blockers encountered during scoped cfd-python verification.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-python`.
- [x] Run `cargo nextest run -p cfd-python --no-tests pass` (0 test binaries).
- [x] Run `cargo fmt -p cfd-python --check`.
- [x] Run cfd-python source/manifest scans for removed num-traits references.
- [x] Run `cargo check -p cfd-1d` and file-level rustfmt for the dependency clippy cleanup.
- [x] Run `cargo clippy -p cfd-python --all-targets -- -D warnings`.
- [x] Run touched-file rustfmt checks and `git diff --check`.

---

# Sprint 1.96.27 Checklist: cfd-io Eunomia Scalar Boundary
**Goal**: Remove cfd-io's direct num-traits dependency by routing checkpoint, binary, and CSV scalar bounds/conversions through Eunomia.

**Success Criteria**:
- ✅ `cfd-io` depends on `eunomia.workspace = true`, not `num-traits.workspace = true`.
- ✅ Checkpoint validation, binary helpers, and CSV helpers use Eunomia `RealField`.
- ✅ `crates/cfd-io` has no `num_traits`, `num-traits`, `FromPrimitive`, or `ToPrimitive` hits.
- ✅ Checkpoint mass-conservation dimension conversion rejects invalid/excessive mesh dimensions instead of silently substituting unit spacing.
- ✅ Existing checkpoint roundtrip/property tests remain value-semantic and green.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct cfd-io `num-traits` use sites in checkpoint, binary, CSV, and manifest.
- [x] Confirm Eunomia exposes `RealField`, `from_f64`, `ZERO`, and finite/absolute-value operations needed by cfd-io.
- [x] Select cfd-io scalar bounds/conversions as the bounded Eunomia migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace `num_traits::Float` bounds with `eunomia::RealField`.
- [x] Replace `FromPrimitive`/`ToPrimitive` conversions in checkpoint mass conservation with Eunomia `from_f64` plus exact `usize` to `f64` validation.
- [x] Remove direct cfd-io `num-traits` dependency and add direct cfd-io `eunomia` dependency.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-io`.
- [x] Run `cargo check -p cfd-io --all-features`.
- [x] Run `cargo nextest run -p cfd-io --no-fail-fast` (3/3 passed).
- [x] Run `cargo clippy -p cfd-io --all-targets --all-features -- -D warnings`.
- [x] Run `cargo fmt -p cfd-io --check`.
- [x] Run cfd-io source/manifest scans for removed num-traits references.

---

# Sprint 1.96.26 Checklist: cfd-io Leto Checkpoint/Binary Boundary
**Goal**: Remove cfd-io's direct and normal-transitive nalgebra dependency path by moving checkpoint and binary dense payloads to Leto arrays and decoupling cfd-io from cfd-core/cfd-math.

**Success Criteria**:
- ✅ `cfd-io` checkpoint fields use Leto `Array2`.
- ✅ Binary vector/matrix helpers use Leto `Array1`/`Array2`.
- ✅ `crates/cfd-io` has no source/manifest `nalgebra`, `DMatrix`, `DVector`, or `RealField` hits.
- ✅ The default normal cfd-io dependency graph has no `nalgebra` package.
- ✅ Checkpoint roundtrip tests preserve shape, row-major values, and checksum equality.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit `cfd-io` direct nalgebra use sites in checkpoint, binary, CSV, tests, and manifest.
- [x] Confirm the pinned Leto revision exposes `Array1`/`Array2`, `shape`, `size`, `get`, `from_elem`, `from_shape_vec`, and `into_vec`; avoid relying on newer local-only iterator APIs.
- [x] Select cfd-io checkpoint/binary dense payloads plus dependency-graph decoupling as the bounded migration slice.

### Phase 2: Execution (10-50%)
- [x] Replace checkpoint `DMatrix` fields with Leto `Array2`.
- [x] Add explicit row-major checkpoint serde payloads because the pinned Leto array type does not derive serde.
- [x] Replace binary `DVector`/`DMatrix` helpers with Leto `Array1`/`Array2`.
- [x] Replace cfd-io CSV `RealField` bounds with `num_traits::Float` as the temporary scalar boundary until Eunomia migration.
- [x] Add a local cfd-io file-format error type and remove unused `cfd-core`/`cfd-math` dependencies.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-io`.
- [x] Run `cargo nextest run -p cfd-io --no-fail-fast` (3/3 passed).
- [x] Run `cargo check -p cfd-io --all-features`.
- [x] Run `cargo fmt -p cfd-io --check`.
- [x] Run `cargo tree -p cfd-io -e normal -i nalgebra` and source scans.
- [x] Record residual Eunomia scalar and optional VTK transitive nalgebra risks.

---

# Sprint 1.96.25 Checklist: cfd-python Leto PyO3 Boundary
**Goal**: Move the bounded `cfd-python` 2D NumPy-return helpers away from direct `ndarray`/`nalgebra` ownership while preserving the Python-visible `numpy.ndarray` return contract.

**Success Criteria**:
- ✅ `cfd-python` no longer depends directly on `ndarray` or `nalgebra`.
- ✅ 2D Poiseuille analytical velocity and Ghia benchmark helper data are constructed as Leto `Array2` values in Rust.
- ✅ NumPy remains only the PyO3 boundary copy target for Python callers.
- ✅ Bounded Cargo check and source scans verify the changed binding path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit direct `ndarray`/`nalgebra` references in `crates/cfd-python`.
- [x] Confirm NumPy's `PyArray2::from_vec2_bound` and Leto's `Array2::from_shape_fn` are available in the locked dependency versions.
- [x] Select the 2D PyO3 array-return helpers as the bounded migration slice.

### Phase 2: Execution (10-50%)
- [x] Add Leto as the cfd-python dense-array dependency.
- [x] Add a private `solver_2d` helper that copies Leto `Array2` rows into NumPy at the PyO3 boundary.
- [x] Replace `DMatrix` and direct `ndarray::Array2` construction in Poiseuille and cavity helper paths.
- [x] Remove direct `ndarray` and `nalgebra` dependencies from `crates/cfd-python/Cargo.toml`.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-python`.
- [x] Run `cargo nextest run -p cfd-python --no-tests pass` (0 test binaries).
- [x] Run `cargo fmt -p cfd-python --check`.
- [x] Run static scans for `ndarray`, `nalgebra`, and `DMatrix` in `crates/cfd-python`.
- [x] Record residual cfd-python `num-traits`/Eunomia risk.

---

# Sprint 1.96.24 Checklist: CFDrs Moirai GPU Boundary
**Goal**: Move the next bounded CFDrs GPU synchronization boundary toward Atlas-owned concurrency while preserving the current WGPU kernel API until the Hephaestus version normalization is complete.

**Success Criteria**:
- ✅ `cfd-core` GPU context creation, GPU support detection, and Poisson residual readback use `moirai::block_on`.
- ✅ `cfd-math` GPU operator sync dispatch and dispatch-metrics test use `moirai::block_on`.
- ✅ The workspace no longer has a `pollster` dependency or source-level `pollster` references.
- ✅ The Hephaestus replacement blocker is recorded with the concrete WGPU version mismatch.
- ✅ Bounded Cargo check and nextest library verification pass for `cfd-core`.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit CFDrs GPU provider usage and direct `pollster` call sites.
- [x] Confirm Hephaestus' WGPU provider uses `wgpu 26.0` while CFDrs still uses `wgpu 0.19`.
- [x] Select the bounded Moirai boundary migration instead of wiring incompatible raw WGPU types.

### Phase 2: Execution (10-50%)
- [x] Replace `pollster::block_on` with `moirai::block_on` in `cfd-core` GPU creation, detection, and residual readback.
- [x] Replace `pollster::block_on` with `moirai::block_on` in `cfd-math` GPU operator sync dispatch and dispatch-metrics test.
- [x] Keep the adapter fallback probe fully async instead of nesting synchronous blockers inside `GpuContext::create_async`.
- [x] Remove `pollster` from `crates/cfd-core/Cargo.toml`, `crates/cfd-math/Cargo.toml`, and the workspace dependency graph.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-core --features gpu`.
- [x] Run `cargo nextest run -p cfd-core --features gpu --lib` (151/151 passed).
- [x] Run `cargo check -p cfd-math --features gpu`.
- [x] Run `cargo nextest run -p cfd-math --features gpu --lib` (279/279 passed).
- [x] Run static scans/checks for workspace `pollster`, touched-file rustfmt, and touched-file diff whitespace.
- [x] Record residual Hephaestus WGPU normalization risk.

---

# Sprint 1.96.23 Checklist: cfd-3d Atlas Array Migration
**Goal**: Move the next bounded `cfd-3d` spectral/NUFFT dense-array path toward Atlas-owned providers without depending on unpublished Apollo worktree state.

**Success Criteria**:
- ✅ `cfd-3d` spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT algorithm paths use Leto arrays for dense storage.
- ✅ Current Apollo FFT/NUFFT `ndarray` API usage is isolated in one private CFDrs boundary module.
- ✅ Bounded Cargo check and nextest library verification pass for `cfd-3d`.
- ✅ Residual dependency on Apollo's published ndarray API is recorded as boundary-only migration risk.

### Phase 1: Foundation & Specs (0-10%)
- [x] Audit CFDrs spectral/NUFFT array ownership and Apollo FFT/NUFFT API shape.
- [x] Confirm the Leto-backed Apollo commit is not reachable from the remote and avoid path-depending on dirty local Apollo state.
- [x] Select the bounded migration slice around `cfd-3d` spectral and IBM NUFFT dense arrays.

### Phase 2: Execution (10-50%)
- [x] Add Leto as a no-default workspace dependency with `ndarray-compat`.
- [x] Add the private `cfd-3d` Apollo array boundary adapter.
- [x] Migrate spectral DNS, forcing, diagnostics, Fourier wrapper, and IBM NUFFT call sites to Leto arrays and checked accessor helpers.

### Phase 3: Closure (50%+)
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo nextest run -p cfd-3d --no-default-features --lib` (206/206 passed).
- [x] Record residual Apollo ndarray boundary risk.

---

# Sprint 1.96.22 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and crates used by `cfd-1d`, then correct the next highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and external graph/algebra/concurrency/serialization crates used by `cfd-1d`.
- ✅ Flow analysis classifies regimes from Reynolds number magnitude under reverse flow.
- ✅ Forward and reverse flow have identical scalar diagnostics where orientation should not affect the scalar quantity.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` flow-analysis path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-1d` flow-analysis code, network graph usage, and direct dependency roles.
- [x] Classify dependency roles for core fluid properties, math/scalar bounds, schematic-derived geometry, graph topology, serialization, and parallel support.
- [x] Identify signed Reynolds classification in reverse-flow analysis as the next highest-risk diagnostics physics gap.

### Phase 2: Execution (10-50%)
- [x] Compute regime-classification Reynolds number from velocity magnitude.
- [x] Add real-network reverse-flow transitional classification test.
- [x] Add forward/reverse scalar diagnostic reciprocity test.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched flow-analysis path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features solver::analysis::analyzers::flow --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features solver::analysis::analyzers::flow` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.21 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and crates used by `cfd-1d`, then correct the next highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and external graph/algebra/concurrency/serialization crates used by `cfd-1d`.
- ✅ Droplet-regime dimensionless groups reject invalid dimensional inputs instead of clamping denominators.
- ✅ Flow-regime classification rejects nonfinite or negative capillary numbers.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` droplet-regime path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-1d` manifest, two-phase physics module, and direct dependency roles.
- [x] Classify dependency roles for core physics errors, math kernels, schematic geometry inputs, graph topology, sparse algebra, serialization, and parallel support.
- [x] Identify denominator clamps in capillary, Weber, and Ohnesorge groups as the next highest-risk two-phase physics gap.

### Phase 2: Execution (10-50%)
- [x] Replace surface-tension denominator clamps with positive finite validation.
- [x] Replace length and density denominator clamps with physical-domain validation.
- [x] Return typed physics errors from droplet-regime dimensionless-group helpers and regime analysis.
- [x] Add value-semantic invalid-input tests.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched droplet-regime path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::droplet_regime --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::droplet_regime` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.20 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and crates used by `cfd-1d`, then correct the next highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and external graph/algebra/concurrency/serialization crates used by `cfd-1d`.
- ✅ Junction-loss resistance derives rheology shear from velocity magnitude rather than signed velocity.
- ✅ Junction-loss resistance uses explicit nonnegative wall shear rate when supplied.
- ✅ Junction-loss resistance rejects negative explicit wall shear rate before non-Newtonian rheology evaluation.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` junction-loss path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-1d` manifest, resistance module graph, and direct dependency roles.
- [x] Classify dependency roles for core rheology/errors, math kernels, schematic topology, graph topology, sparse algebra, serialization, and parallel support.
- [x] Identify signed junction-loss shear rate and ignored explicit shear as the next highest-risk resistance-physics gap.

### Phase 2: Execution (10-50%)
- [x] Derive default junction wall shear rate from velocity magnitude.
- [x] Route explicit nonnegative shear rate into junction apparent-viscosity evaluation.
- [x] Reject negative explicit wall shear rate as a physical invariant violation.
- [x] Add value-semantic Casson blood tests for shear-thinning dependence and reverse-flow reciprocity.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched junction-loss resistance path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::junction_loss --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::junction_loss` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.19 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and crates used by `cfd-1d`, then correct the next highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and external graph/algebra/concurrency/serialization crates used by `cfd-1d`.
- ✅ Darcy-Weisbach resistance uses explicit nonnegative wall shear rate when supplied.
- ✅ Darcy-Weisbach resistance rejects negative explicit wall shear rate before non-Newtonian rheology evaluation.
- ✅ Darcy-Weisbach resistance propagates rheology errors rather than falling back to baseline viscosity.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` Darcy-Weisbach path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-1d` manifest, resistance module graph, and direct dependency roles.
- [x] Classify dependency roles for core rheology/errors, math kernels, schematic topology, graph topology, sparse algebra, serialization, and parallel support.
- [x] Identify ignored explicit Darcy-Weisbach shear rate and masked rheology errors as the next highest-risk resistance-physics gap.

### Phase 2: Execution (10-50%)
- [x] Route explicit nonnegative shear rate into Darcy-Weisbach apparent-viscosity evaluation.
- [x] Reject negative explicit wall shear rate as a physical invariant violation.
- [x] Replace baseline-viscosity fallback with direct propagation of `cfd-core` rheology errors.
- [x] Add value-semantic Casson blood tests for shear-thinning dependence and reverse-flow reciprocity.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched Darcy-Weisbach resistance path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::darcy_weisbach --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::darcy_weisbach` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.18 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and crates used by `cfd-1d`, then correct the next highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and external graph/algebra/concurrency crates used by `cfd-1d`.
- ✅ Rectangular-channel resistance uses explicit nonnegative wall shear rate when supplied.
- ✅ Rectangular-channel resistance rejects negative explicit wall shear rate before non-Newtonian rheology evaluation.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` rectangular resistance path.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-1d` manifest, public module graph, and direct dependency roles.
- [x] Classify dependency roles for core rheology/errors, math kernels, schematic topology, graph topology, sparse algebra, serialization, and parallel support.
- [x] Identify ignored explicit rectangular-channel shear rate as the next highest-risk local resistance-physics gap.

### Phase 2: Execution (10-50%)
- [x] Route explicit nonnegative shear rate into rectangular-channel apparent-viscosity evaluation.
- [x] Reject negative explicit wall shear rate as a physical invariant violation.
- [x] Add value-semantic Casson blood tests for shear-thinning dependence and reverse-flow reciprocity.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched rectangular resistance path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::rectangular --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::rectangular` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.17 Checklist: cfd-3d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-3d` and direct dependencies, then correct the next highest-risk cfd-3d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-3d`, `cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates.
- ✅ FEM Stokes validation rejects invalid fluid density and viscosity before assembly.
- ✅ FEM Stokes validation rejects invalid pressure-space, body-force, boundary-data, and element-viscosity states.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-3d` FEM paths.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-3d` manifest, public module graph, and direct dependency roles.
- [x] Classify dependency roles for core fluid/BC contracts, math kernels, mesh authority, I/O, lower-fidelity references, schematic topology, and numerical support crates.
- [x] Identify incomplete `StokesFlowProblem::validate` physical invariant checks as the next highest-risk local FEM gap.

### Phase 2: Execution (10-50%)
- [x] Add positive finite density and viscosity validation.
- [x] Add pressure corner-node count validation.
- [x] Add finite body-force and boundary-condition data validation.
- [x] Add per-element viscosity field length and positivity validation.
- [x] Add value-semantic FEM problem validation tests.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched FEM problem path.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features fem::problem --lib`.
- [x] Run `cargo nextest run -p cfd-3d --lib --no-default-features fem::problem` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.16 Checklist: cfd-3d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-3d` and direct dependencies, then correct the next highest-risk cfd-3d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-3d`, `cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates.
- ✅ Level-set Hamilton-Jacobi transport rejects nonpositive and nonfinite time steps.
- ✅ Level-set transport rejects nonpositive/nonfinite grid spacing and nonfinite velocity components before derivative reconstruction.
- ✅ Bounded Cargo check, integration test, clippy, and nextest verification pass for the touched `cfd-3d` level-set paths.

### Phase 1: Foundation & Specs (0-10%)
- [x] Re-audit `cfd-3d` manifest, public module graph, and direct dependency roles.
- [x] Classify dependency roles for core physics/errors, math kernels, mesh authority, I/O, lower-fidelity references, schematic topology, and numerical support crates.
- [x] Identify missing level-set Hamilton-Jacobi transport precondition checks as the next highest-risk local physics gap.

### Phase 2: Execution (10-50%)
- [x] Reject nonpositive and nonfinite level-set time steps before transport.
- [x] Reject nonpositive/nonfinite grid spacing before transport.
- [x] Reject nonfinite velocity components before WENO or first-order derivative reconstruction.
- [x] Add value-semantic invalid-input tests for level-set transport.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched level-set paths.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features --test level_set_tests`.
- [x] Run `cargo nextest run -p cfd-3d --test level_set_tests --no-default-features` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.15 Checklist: cfd-3d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-3d` and direct dependencies, then correct the highest-risk cfd-3d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-3d`, `cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, `cfd-2d`, `cfd-schematics`, and numerical support crates.
- ✅ VOF explicit interface transport rejects nonpositive and nonfinite time steps.
- ✅ VOF CFL evaluation rejects nonfinite velocity components before flux calculation.
- ✅ Bounded Cargo check, integration test, clippy, and nextest verification pass for the touched `cfd-3d` VOF paths.

### Phase 1: Foundation & Specs (0-10%)
- [x] Map `cfd-3d` manifest, public module graph, and direct dependencies.
- [x] Classify dependency roles: core physics/errors, math kernels, mesh authority, I/O adapters, lower-fidelity references, schematic topology, FFT/array/concurrency support.
- [x] Identify missing VOF explicit-transport precondition checks as the highest-risk local physics gap.

### Phase 2: Execution (10-50%)
- [x] Reject nonpositive and nonfinite VOF time steps before CFL evaluation.
- [x] Reject nonfinite VOF velocity components before geometric or algebraic advection.
- [x] Add value-semantic invalid-input tests for VOF transport.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched VOF paths.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features --test vof_tests`.
- [x] Run `cargo nextest run -p cfd-3d --test vof_tests --no-default-features` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.14 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and direct dependencies, then correct the highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and how each dependency is used.
- ✅ Hagen-Poiseuille shear-rate derivation is invariant under flow reversal for non-Newtonian fluids.
- ✅ Explicitly negative shear-rate inputs are rejected instead of being passed to rheology models.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` module.

### Phase 1: Foundation & Specs (0-10%)
- [x] Map `cfd-1d` manifest, public module graph, and direct internal dependencies.
- [x] Classify dependency roles: core rheology/errors, numerical kernels, schematic topology and cross-section input.
- [x] Identify signed derived shear rate in Hagen-Poiseuille as the highest-risk local physics gap.

### Phase 2: Execution (10-50%)
- [x] Use velocity magnitude when deriving wall shear rate from velocity or flow rate.
- [x] Reject explicit negative wall shear-rate inputs.
- [x] Add Casson blood reverse-flow reciprocity and negative-shear tests.

### Phase 3: Closure (50%+)
- [x] Run marker scan for the touched Hagen-Poiseuille path.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::hagen_poiseuille --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::hagen_poiseuille` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.13 Checklist: cfd-2d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-2d` and direct dependencies, then correct the highest-risk cfd-2d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-2d`, `cfd-core`, `cfd-math`, `cfd-mesh`, `cfd-io`, `cfd-1d`, and `cfd-schematics`.
- ✅ LBM initialization rejects velocities outside the low-Mach incompressible regime.
- ✅ LBM Zou-He velocity and pressure boundary reconstruction reject derived or prescribed `Ma > 0.1`.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-2d` LBM paths.

### Phase 1: Foundation & Specs (0-10%)
- [x] Map `cfd-2d` manifest, public module graph, and direct internal dependencies.
- [x] Classify dependency roles: core physics/errors, numerical kernels, mesh and IO adapters, 1D reference seeding, schematic topology.
- [x] Identify missing low-Mach enforcement in D2Q9 LBM velocity entry points.

### Phase 2: Execution (10-50%)
- [x] Add low-Mach validation for LBM initialization velocities.
- [x] Add low-Mach validation for Zou-He velocity boundaries and pressure-derived velocities.
- [x] Return boundary validation errors through `LbmSolver::step`.
- [x] Add value-semantic high-Mach rejection tests.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched cfd-2d LBM paths.
- [x] Run `cargo check -p cfd-2d --no-default-features`.
- [x] Run `cargo test -p cfd-2d --no-default-features solvers::lbm --lib`.
- [x] Run `cargo nextest run -p cfd-2d --lib --no-default-features solvers::lbm` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-2d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.12 Checklist: cfd-1d Dependency-Aware Physics Audit
**Goal**: Audit `cfd-1d` and direct dependencies, then correct the highest-risk cfd-1d physics gap found.

**Success Criteria**:
- ✅ Audit covers `cfd-1d`, `cfd-core`, `cfd-math`, `cfd-schematics`, and how each dependency is used.
- ✅ Serpentine scalar resistance losses are invariant under flow reversal.
- ✅ Tests inspect coefficients, resistance, Reynolds number, wall shear rate, Dean number, and pressure-drop magnitudes.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched `cfd-1d` module.

### Phase 1: Foundation & Specs (0-10%)
- [x] Map `cfd-1d` manifest, public module graph, and direct internal dependencies.
- [x] Classify dependency roles: core physics/errors, numerical linear algebra, schematic topology authority.
- [x] Identify reverse-flow sign leakage in serpentine scalar-loss physics.

### Phase 2: Execution (10-50%)
- [x] Use `|u|` for serpentine shear-rate, Reynolds, Dean, friction, and bend-loss calculations.
- [x] Preserve resistance-coefficient invariance under reversed velocity.
- [x] Add value-semantic reverse-flow tests.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched cfd-1d paths.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::serpentine --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::serpentine` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-1d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.11 Checklist: cfd-3d Sigma SGS Energy Physics
**Goal**: Correct `cfd-3d` Sigma LES turbulent kinetic-energy diagnostics.

**Success Criteria**:
- ✅ Sigma SGS kinetic energy uses `k_sgs = (nu_t / (C_k Delta))^2`.
- ✅ Sigma shares the same SGS energy conversion used by WALE and Vreman.
- ✅ Tests inspect the computed center-cell energy and distinguish the former strain-rate formula.
- ✅ Bounded Cargo check, unit test, clippy, and nextest verification pass for the touched Sigma module.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify dimensional mismatch in `SigmaModel::turbulent_kinetic_energy`.
- [x] Specify the Yoshizawa invariant: eddy viscosity `nu_t = C_k Delta sqrt(k_sgs)`.

### Phase 2: Execution (10-50%)
- [x] Replace Sigma-specific `nu_t |S| / Delta` energy computation.
- [x] Route Sigma through shared `kinetic_energy_from_eddy_viscosity`.
- [x] Add value-semantic regression coverage for the corrected relation.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched Sigma paths.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features physics::turbulence::sigma --lib`.
- [x] Run `cargo nextest run -p cfd-3d --lib --no-default-features physics::turbulence::sigma` under a 30-second shell timeout.
- [x] Run `cargo clippy -p cfd-3d --no-default-features --lib -- -W clippy::all -W clippy::pedantic`.

---

# Sprint 1.96.10 Checklist: cfd-3d VOF Directional CFL Physics
**Goal**: Correct `cfd-3d` VOF timestep selection for diagonal and anisotropic flow.

**Success Criteria**:
- ✅ `VofSolver::calculate_timestep` uses the summed directional advective rate enforced by VOF advection.
- ✅ Diagonal flow with target CFL 1 returns `dt = 1/3` on a unit grid.
- ✅ The former norm/min-spacing estimate is rejected by the geometric VOF CFL check.
- ✅ Bounded Cargo check, integration test, and nextest verification pass for the touched VOF target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify mismatch between solver timestep selection and `AdvectionMethod` CFL enforcement.
- [x] Specify the invariant: `max_cells(|u_x|dt/dx + |u_y|dt/dy + |u_z|dt/dz) <= C`.

### Phase 2: Execution (10-50%)
- [x] Replace Euclidean speed/min-spacing timestep with reciprocal summed directional advective rate.
- [x] Add Rustdoc theorem and proof sketch to the timestep API.
- [x] Add value-semantic diagonal-flow regression coverage.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched VOF paths.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features --test vof_tests test_calculate_timestep_uses_directional_vof_cfl`.
- [x] Run `cargo nextest run -p cfd-3d --test vof_tests --no-default-features test_calculate_timestep_uses_directional_vof_cfl` under a 30-second shell timeout.

---

# Sprint 1.96.9 Checklist: cfd-2d Upwind Coefficient Orientation
**Goal**: Correct `cfd-2d` first-order upwind finite-volume coefficients for west-face flow orientation.

**Success Criteria**:
- ✅ East-face coefficient remains `a_E = D_E + max(-F_E, 0)`.
- ✅ West-face coefficient uses `a_W = D_W + max(F_W, 0)`.
- ✅ Tests inspect positive and negative west-face flux values.
- ✅ Bounded Cargo check, unit test, and nextest verification pass for the touched library target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify sign mismatch in `FirstOrderUpwind::coefficients` for west-face convection.
- [x] Specify the coefficient invariant from finite-volume upwinding: neighbor coefficients must remain nonnegative and select the upstream cell by face-flow orientation.

### Phase 2: Execution (10-50%)
- [x] Correct `a_W` to use `max(F_W, 0)`.
- [x] Update the upwind theorem docs with east/west coefficient form.
- [x] Add value-semantic tests for west-face positive and negative flow.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched convection paths.
- [x] Run `cargo check -p cfd-2d --no-default-features`.
- [x] Run `cargo test -p cfd-2d --no-default-features discretization::convection --lib`.
- [x] Run `cargo nextest run -p cfd-2d --lib --no-default-features discretization::convection` under a 30-second shell timeout.

---

# Sprint 1.96.8 Checklist: cfd-1d Branch Reverse-Flow Physics
**Goal**: Correct `cfd-1d` branch-junction solvers for reversed parent flow orientation.

**Success Criteria**:
- ✅ Two-way and three-way pressure-balanced branch solvers accept negative parent flow.
- ✅ Prescribed split solvers preserve negative daughter-flow orientation.
- ✅ Wall-shear and apparent viscosity remain nonnegative magnitude diagnostics.
- ✅ Pressure drops reverse sign with flow orientation.
- ✅ Bounded Cargo check, unit test, and nextest verification pass for the touched library target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify sign-sensitive branch-junction checks that rejected `Q_parent < 0`.
- [x] Specify the orientation invariant: reversing all branch flows changes signed flows and pressure drops but leaves shear-rate and viscosity magnitudes unchanged.

### Phase 2: Execution (10-50%)
- [x] Solve pressure-balanced splits on `|Q_parent|` and reapply parent-flow orientation to daughter flows.
- [x] Remove nonnegative-flow rejection from prescribed split paths.
- [x] Use `|Q|` for wall-shear and apparent-viscosity inputs.
- [x] Preserve signed pressure-drop calculation through signed `Q`.
- [x] Add value-semantic regression tests for two-way and three-way reversed flow.

### Phase 3: Closure (50%+)
- [x] Run marker scan for touched branch-junction paths.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features --test branch_reverse_flow_orientation`.
- [x] Run `cargo nextest run -p cfd-1d --test branch_reverse_flow_orientation --no-default-features` under a 30-second shell timeout.

---

# Sprint 1.96.7 Checklist: cfd-3d Venturi Pressure-Coefficient Physics
**Goal**: Correct `cfd-3d` Venturi pressure coefficients to use throat dynamic-pressure scaling.

**Success Criteria**:
- ✅ `cp_throat` and `cp_recovery` use `0.5 ρ (Q/A_throat)^2` as the denominator.
- ✅ Coefficient scaling rejects non-positive throat flow or area instead of producing dimensionless values from an undefined scale.
- ✅ New tests validate computed coefficient values and the zero-flux rejection path.
- ✅ Bounded Cargo check, unit test, and nextest verification pass for the touched library target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify mismatch between the Venturi pressure-recovery theorem and the solver's inlet dynamic-pressure denominator.
- [x] Specify the coefficient invariant: pressure coefficients must be nondimensionalized by throat dynamic pressure, not inlet dynamic pressure.

### Phase 2: Execution (10-50%)
- [x] Add a single coefficient helper derived from face-integrated throat flux and throat area.
- [x] Route `cp_throat` and `cp_recovery` through the throat dynamic-pressure helper.
- [x] Document solution coefficient fields as throat-scaled quantities.
- [x] Add value-semantic tests for coefficient values and rejection of undefined dynamic pressure.

### Phase 3: Closure (50%+)
- [x] Run marker scan for coefficient implementation paths.
- [x] Run `cargo check -p cfd-3d --no-default-features`.
- [x] Run `cargo test -p cfd-3d --no-default-features venturi_pressure_coefficients --lib`.
- [x] Run `cargo nextest run -p cfd-3d --lib --no-default-features venturi_pressure_coefficients` under a bounded shell timeout.

---

# Sprint 1.96.6 Checklist: cfd-2d Explicit Stability Physics
**Goal**: Correct `cfd-2d` explicit time-step bounds for 2D advection-diffusion.

**Success Criteria**:
- ✅ `max_stable_dt` uses the summed 2D advection CFL rate `|u|/dx + |v|/dy`.
- ✅ `max_stable_dt` uses the summed 2D diffusion rate `ν(1/dx² + 1/dy²)`.
- ✅ New tests validate the returned `dt` against `advection_cfl` and `diffusion_number`.
- ✅ Bounded Cargo check, unit test, and nextest verification pass for the touched library target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify that componentwise `min(dx/|u|, dy/|v|)` can allow summed 2D CFL greater than 1.
- [x] Identify that `0.5*min(dx²,dy²)/ν` overestimates the 2D explicit diffusion bound on square grids by a factor of 2.

### Phase 2: Execution (10-50%)
- [x] Implement reciprocal summed-rate advection and diffusion time-step bounds.
- [x] Document the stability theorem and proof sketch in the implementation.
- [x] Add value-semantic tests for anisotropic-grid advection and diffusion cases.

### Phase 3: Closure (50%+)
- [x] Run marker scan for approximation/stub wording in the touched CFL module.
- [x] Run `cargo check -p cfd-2d --no-default-features`.
- [x] Run `cargo test -p cfd-2d --no-default-features stability::cfl --lib`.
- [x] Run `cargo nextest run -p cfd-2d --lib --no-default-features stability::cfl` under a 30-second shell timeout.

---

# Sprint 1.96.5 Checklist: cfd-1d Venturi Reverse-Flow Physics
**Goal**: Correct `cfd-1d` Venturi resistance and analysis for reverse-flow inputs.

**Success Criteria**:
- ✅ Negative inlet velocity no longer routes coefficient calculation through the zero-flow branch.
- ✅ Reynolds number, shear rate, viscosity query, friction factor, and scalar pressure-loss terms use velocity magnitude.
- ✅ Reported throat velocity preserves the input flow orientation.
- ✅ Symmetric Venturi scalar resistance coefficients are invariant under velocity sign reversal.
- ✅ Bounded Cargo check, unit test, and nextest verification pass for the touched library target.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify sign-sensitive Venturi branches in `calculate_coefficients` and `analyze`.
- [x] Specify the symmetric-geometry invariant: reversing flow orientation changes signed velocities but not scalar loss magnitudes.

### Phase 2: Execution (10-50%)
- [x] Add a local magnitude helper for generic scalar values.
- [x] Use `|V_inlet|`, `|V_throat|`, and `|Q|` for coefficient decomposition.
- [x] Use shear-rate and Reynolds magnitudes in detailed Venturi analysis.
- [x] Add reverse-flow regression tests for coefficients, resistance, and pressure decomposition.

### Phase 3: Closure (50%+)
- [x] Run marker scan for approximation/stub wording in the touched venturi directory.
- [x] Run `cargo check -p cfd-1d --no-default-features`.
- [x] Run `cargo test -p cfd-1d --no-default-features physics::resistance::models::venturi --lib`.
- [x] Run `cargo nextest run -p cfd-1d --lib --no-default-features physics::resistance::models::venturi` under a 30-second shell timeout.

---

# Sprint 1.96.4 Checklist: Geometric Conservation Residual Verification
**Goal**: Replace copy-through geometric conservation checks with residual-based Euler and Runge-Kutta updates.

**Success Criteria**:
- ✅ `cfd-validation` GCL checks evaluate a conservative second-order finite-volume residual.
- ✅ Constant fields are preserved by Euler, midpoint, SSPRK3, and RK4 stage formulas.
- ✅ Non-constant quadratic fields produce the analytically expected residual and state update.
- ✅ Unsupported RK stage counts return a typed rejection.
- ✅ Bounded Cargo check, unit test, and nextest verification pass within 30 seconds after target scoping.

### Phase 1: Foundation & Specs (0-10%)
- [x] Identify copy-through Euler/RK evolution in `crates/cfd-validation/src/conservation/geometric.rs`.
- [x] Specify the constant-state invariant: conservative face-flux gradients vanish exactly for `u_ij = c`.

### Phase 2: Execution (10-50%)
- [x] Add `conservative_residual` using east/west/north/south face-flux divergence.
- [x] Replace Euler copy-through with `u^{n+1} = u^n + dt R(u^n)`.
- [x] Replace RK copy-through with residual-based midpoint, SSPRK3, and RK4 stage formulas.
- [x] Add value-semantic regression tests for quadratic residual response and unsupported RK stages.

### Phase 3: Closure (50%+)
- [x] Run marker scan for copy-through and simplification wording in the touched GCL module.
- [x] Run `cargo check -p cfd-validation --no-default-features`.
- [x] Run `cargo test -p cfd-validation --no-default-features conservation::geometric --lib`.
- [x] Run `cargo nextest run -p cfd-validation --lib --no-default-features conservation::geometric` under a 30-second shell timeout.

---

# Sprint 1.96.3 Checklist: Womersley Analytical SSOT
**Goal**: Replace validation-local Womersley approximations with the canonical exact complex-Bessel implementation and verify no-slip behavior.

**Success Criteria**:
- ✅ `cfd-validation` Womersley velocity, wall shear stress, and flow rate use `cfd-1d` `WomersleyProfile`.
- ✅ Rustdoc documents the exact Womersley Bessel solution and no-slip proof sketch.
- ✅ Wall no-slip regression checks computed velocity values at multiple phases.
- ✅ Bounded check and unit-test verification pass.
- ⚠️ Targeted nextest exceeded the 60-second compile bound before test execution.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `cfd-validation` depends on `cfd-1d`, so reusing the canonical Womersley profile introduces no dependency cycle.
- [x] Identify validation-local approximate velocity, wall-stress, and flow-rate formulas.

### Phase 2: Execution (10-50%)
- [x] Add a private exact-profile adapter for the validation `WomersleyFlow`.
- [x] Delegate velocity, wall shear stress, and flow rate to `cfd-1d` exact Womersley profile.
- [x] Add a no-slip wall velocity regression over multiple phases.

### Phase 3: Closure (50%+)
- [x] Run marker scan for approximate/simplified wording in the touched Womersley module.
- [x] Run `cargo check -p cfd-validation --no-default-features`.
- [x] Run `cargo test -p cfd-validation --no-default-features analytical::womersley --lib`.
- [x] Record targeted nextest compile-bound timeout.

---

# Sprint 1.96.2 Checklist: Optimization Terminology Contract Cleanup
**Goal**: Remove misleading unsupported/stub-like wording from boundary and serpentine optimization contracts without changing numerical behavior.

**Success Criteria**:
- ✅ `cfd-schematics` serpentine optimization generator is named for its role in objective evaluation.
- ✅ Boundary stencil rejection message states unsupported order instead of incomplete implementation.
- ✅ Stale symbol and marker scan for touched paths is clean.
- ✅ Bounded Cargo check passes for touched crates.
- ⚠️ Targeted nextest exceeded the 60-second compile bound before test execution.

### Phase 1: Foundation & Specs (0-10%)
- [x] Confirm `generate_simplified_serpentine_path` is `pub(super)` internal to `cfd-schematics` geometry optimization.
- [x] Confirm boundary unsupported-order change affects diagnostic wording only.

### Phase 2: Execution (10-50%)
- [x] Rename `generate_simplified_serpentine_path` to `generate_optimization_serpentine_path`.
- [x] Update Nelder-Mead and grid-search call sites.
- [x] Replace `"Stencil order ... not implemented"` with `"Stencil order ... is unsupported"`.

### Phase 3: Closure (50%+)
- [x] Run stale-symbol scan for renamed function and touched marker terms.
- [x] Run `cargo check -p cfd-core -p cfd-schematics --no-default-features`.
- [x] Record targeted nextest compile-bound timeout.

---

# Sprint 1.96.1 Checklist: Workspace SSOT Cleanup
**Goal**: Remove obsolete tracked root-source artifacts that duplicate canonical crate implementations and are not reachable from Cargo, tests, examples, docs, or report assets.

**Success Criteria**:
- ✅ Unreferenced root `old_*.rs` historical source files and empty `csg_bi.rs` are removed.
- ✅ Repository reference scan confirms no authoritative artifact imports or cites the removed files.
- ✅ Bounded Cargo verification completes for the affected root package metadata path.
- ✅ Windows GNU builds use MSYS2 clang and lld instead of gcc.
- ✅ Misleading marker terminology is removed from selected explicit unsupported-operation paths and bounded-model comments.
- ⚠️ Workspace `cargo nextest` was attempted after clang/lld configuration and exceeded the 60-second compile bound before test execution.

### Phase 1: Foundation & Specs (0-10%)
- [x] Classify the cleanup as a patch-class SSOT artifact removal with no public API, algorithm, or report-output change.
- [x] Verify the files are tracked and not referenced by Cargo manifests, crate modules, examples, tests, docs, or report assets.

### Phase 2: Execution (10-50%)
- [x] Delete `old_assemble.rs`, `old_arrangement.rs`, `old_phase2.rs`, `old_operations.rs`, `old_indexed.rs`, `old_gwn_bvh.rs`, `old_seam.rs`, `old_phase4.rs`, and `csg_bi.rs`.
- [x] Configure `x86_64-pc-windows-gnu` Cargo builds to link through MSYS2 `clang.exe` with `-fuse-ld=lld`.
- [x] Configure C/C++ build-script tools to use MSYS2 `clang.exe`, `clang++.exe`, `llvm-ar.exe`, and `llvm-ranlib.exe`.
- [x] Replace misleading `mock`, `placeholder`, `not implemented`, and `simplified` wording in selected code paths where the executable contract is explicit.
- [x] Update backlog and gap-audit artifacts for traceability.

### Phase 3: Closure (50%+)
- [x] Run bounded verification for reference absence and Cargo metadata/build impact.
- [x] Record workspace nextest compile-time timeout separately from cleanup correctness.
- [x] Run bounded package check for touched crates after marker cleanup.
- [x] Record targeted nextest compile-time timeout separately from source correctness.

---

# Sprint 1.96.0 Checklist: HCOC Cellular Injury & CTC Detection
**Goal**: Integrate mathematical models for cellular cavitation-induced injury and stiffness-coupled anomalous nucleation into the core engine.

**Success Criteria**:
- ✅ Mathematical proofs provided for cell failure grading and heterogeneous inception.
- ✅ `bio_damage` module verified via Proptest, demonstrating rigorous threshold invariants.
- ✅ `heterogeneous_nucleation` matches literature inception divergence ratios between cell types.

### Phase 1: Foundation & Specs (0-10%)
- [x] Specify mathematical domain representing membrane shear stress and threshold porosity limits. 

### Phase 2: Execution (10-50%)
- [x] Implement `bio_damage.rs` evaluating cavitation-induced cell injury fractions from Rayleigh collapse loading and ordered membrane strain thresholds.
- [x] Implement `heterogeneous_nucleation.rs` extending nuclei scalar fields using physical Blake bounds.
- [x] Property test implementation against derived constraints.
- [x] Fixed E0282 broken `apollofft` test blocking test runner in CI.

### Phase 3: Closure (50%+)
- [x] Sync documentation rules and examples.
- [x] Couple nuclei diffusion into the 3D cavitation transport solver and verify advective-diffusive spreading.
- [x] Reuse the 3D cavitation-source workspace to remove the per-step source allocation hot path.
- [x] Validate 3D cavitation flow-field dimensions before stepping to prevent panic-only failure paths.
- [x] Validate pressure and density dimensions in 3D cavitation damage accumulation and remove repeated matrix indexing from the hot path.
- [x] Validate cavitation-source dimensions in the 3D cavitation transport helper and clamp source updates to feasible bounds.
- [x] Validate nuclei transport dimensions in the 3D cavitation solver and use slice-based advection-diffusion accumulation.
- [x] Replace schematic auto-layout map churn with an indexed borrowed layout cache and index-keyed parallel channel grouping.
- [x] Close the `cfd-schematics` geometry-bounds review with canonical constant ranges and clearance-width relation validation.
- [x] Unify Milestone 12 GA convergence narrative and figure annotations on the trailing fitness window trend.
- [x] Harden `cfd-ui` clipping commands to preserve slot identity and reject stale undo restores.
- [x] Add Milestone 12 validation evidence and artifact traceability to the report, results artifact, and asset-review manifest.
- [x] Make the Milestone 12 release report emit the final authoritative narrative in one pass so the example completes after review.
- [x] Bind the Apollo-backed periodic DNS stepper to a validated reusable `FftPlan3D`.
- [x] Replace the approximate `cfd-2d` MUSCL3/QUICK face reconstruction with exact quadratic interpolation and regression tests.
- [x] Replace simplified `cfd-1d` margination lift aggregation with separated wall-induced and shear-gradient inertial scaling.
- [x] Make `cfd-1d` droplet occupied-channel snapshots a finite-length occupancy projection with regression coverage.
- [x] Replace the `cfd-2d` turbulence benchmark placeholder branch with typed supported-model dispatch.
- [x] Resolve the coupled pressure-event blood hematocrit regression with finite startup viscosity, row-equilibrated pressure solves, and duplicate-entry-preserving dense fallback conversion.
- [x] Replace compact `cfd-1d` plasma-skimming screening with threshold-aware Pries phase separation using Murray-inferred sibling geometry.
- [x] Replace `cfd-2d` serpentine mixing's exponential estimate with the closed-form transverse diffusion eigenfunction series and report analytical L90/t90 in the discretized solver.
- [x] Replace `cfd-3d` LES turbulent-kinetic-energy viscosity aliases with a shared Yoshizawa SGS energy relation and regression tests.
- [x] Remove silent clamps and caps from `cfd-2d` Pries plasma-skimming by adding a checked physical-envelope evaluator and value-semantic threshold tests.
- [x] Replace `cfd-2d` WALE boundary zero-gradient assumptions with second-order one-sided derivative stencils and polynomial reproduction tests.
- [x] Replace `cfd-3d` Spalart-Allmaras all-zero TKE with a Yoshizawa wall-distance diagnostic and regression tests.
- [x] Reject uninitialized `cfd-3d` k-epsilon state rather than returning synthetic zero viscosity, TKE, or dissipation fields.
- [x] Replace `cfd-2d` Smagorinsky LES zero SGS energy/dissipation placeholders, boundary zero-strain enforcement, and default SGS viscosity floor with documented diagnostics and value-semantic tests.
- [x] Close review finding 1 by replacing the `cfd-1d` margination singular wall-lift cutoff and public clamp with an explicit validated `[0, 1]` envelope and derived force regression tests.
- [x] Confirm review finding 2 remains closed by typed `cfd-2d` turbulence benchmark dispatch with unsupported-model rejection before benchmark execution.
- [x] Confirm review finding 3 remains closed by finite-span-derived droplet occupied-channel projection and consistency tests.
- [x] Close review finding 3 at the representation level by removing stored point-droplet occupied-channel state and deriving occupied channels from finite-length spans.
- [x] Remove residual nonzero Smagorinsky SGS floors from `cfd-2d` turbulence validation configurations.
- [x] Correct Milestone 12 cross-mode therapy utility so Option 1 receives acoustic-cavitation credit from ultrasound resonance instead of being capped as separation-only.

---

# Sprint 1.95.1 Checklist: CFD-MESH 3D Performance Optimization

## Sprint Overview
**Goal**: Resolve WASM OOM and execution panics by flattening memory arrays and ensuring zero-allocation loops in 3D Delaunay generation.

**Success Criteria**:
- ✅ `BowyerWatson3D::insert_point` executes with 0 heap allocations per point.
- ✅ Delaunay sphere WASM generation at `res=0.15` completes successfully without `unreachable` panic.
- ✅ 100% of property tests pass in `cfd-mesh`.
- ✅ No approximations or empirical epsilon tuning.

## Current Sprint Tasks

### Phase 1: Foundation & Audit (0-10%)
- [x] **Task 1.1**: Audit WASM OOM root cause.
  - [x] Identify O(N²) quadratic allocation overhead in `BowyerWatson3D::insert_point`.
  - [x] Create gap analysis for memory structures.

### Phase 2: Execution (10-50%)
- [x] **Task 2.1**: Refactor `BowyerWatson3D` Memory
  - [x] Add `cavity_cache: HashMap` and `next_tetrahedra: Vec` inside the `BowyerWatson3D` struct.
  - [x] Use `clear()` and variable swapping to prevent reallocations.
  - [x] Implement mathematical proof of invariant preservation in docs.
- [x] **Task 2.2**: Optimize `SdfMesher::build_volume`
  - [x] Pre-allocate `bwid_to_vid`, `used`, and `mesh` internals using exact point counting.
  - [x] Test bounding-box culling efficiency.

### Phase 3: Verification (50%+)
- [x] **Task 3.1**: Property Validation
  - [x] Execute `cargo nextest run -p cfd-mesh`.
  - [x] Verify Euler-Poincaré invariant on coarse and fine meshes.
- [x] **Task 3.2**: WASM End-to-End
  - [x] Rebuild `cfd-ui` WASM target.
  - [x] Confirm generation completes in the browser at fine resolutions.

---

# Sprint 1.91.0 Checklist: Advanced Validation Framework Expansion

### Phase 1: MMS Framework Extension (Week 1-2)
- [x] **Task 1.1**: Design MMS geometry abstraction layer ✅ COMPLETED
  - [x] Define geometry interface for MMS source term generation
  - [x] Implement coordinate transformation system
  - [x] Add geometry validation and boundary handling
  - [x] Implement RectangularDomain geometry with boundary conditions
  - [x] Add comprehensive tests and documentation
- [x] **Task 1.2**: Implement circular domain MMS ✅ COMPLETED
  - [x] Create circular geometry class with boundary detection
  - [x] Implement CircularDomain with full Geometry trait support
  - [x] Add boundary normal calculation and parametric coordinates
  - [x] Comprehensive test suite for circular domain operations
- [x] **Task 1.3**: Implement annular domain MMS ✅ COMPLETED
  - [x] Extend circular geometry for annular regions
  - [x] Implement AnnularDomain with full Geometry trait support
  - [x] Handle inner/outer boundary conditions with separate normal calculations
  - [x] Comprehensive test suite for annular domain operations
  - [x] Validate MMS accuracy for annular geometries with proper area calculations

### Phase 2: Richardson Extrapolation (Week 3-4)
- [x] **Task 2.1**: Core Richardson extrapolation library ✅ COMPLETED
  - [x] Implement grid refinement algorithms
  - [x] Add error estimation and convergence rate calculation
  - [x] Create extrapolation result data structures
- [x] **Task 2.2**: Integration with MMS framework ✅ COMPLETED
  - [x] Connect Richardson extrapolation to MMS solvers
  - [x] Implement automated grid convergence studies
  - [x] Add convergence plotting and analysis
- [x] **Task 2.3**: Validation and testing ✅ COMPLETED
  - [x] Test Richardson extrapolation accuracy
  - [x] Validate convergence rate estimation
  - [x] Add comprehensive test suite

### Phase 2b: Architectural Integrity Remediation (Emergency Audit)
- [x] **Task 2.4**: Resolve Critical Audit Gaps ✅ COMPLETED
  - [x] Fix redundant ILU implementations (Removed ILUPreconditioner)
  - [x] Rename deceptive SchwarzPreconditioner to SerialSchwarzPreconditioner
  - [x] Remove unsubstantiated parallel scalability claims in CG solver
  - [x] Fix flawed DeflationPreconditioner test case
  - [x] **Task 2.4b**: Resolve CFD-CORE Audit Gaps ✅ COMPLETED
    - [x] Fix Fake Distributed GMRES (Implemented Givens/Least Squares)
    - [x] Fix Fake Additive Schwarz (Implemented Local Matrix Assembly/LU)
    - [x] Fix Fake Block Jacobi (Implemented Diagonal Extraction)
  - [x] **Task 2.4c**: Resolve CFD-MESH Audit Gaps ✅ COMPLETED
    - [x] Fix Fake Mesh Refinement (Marked as NotImplemented)
    - [x] Fix Missing Distributed Mesh Support (Added global_id/partition_id)

### Phase 2c: Numerical Correctness Remediation (Legacy Tests)
- [x] **Task 2.5**: Fix Legacy Test Failures ✅ COMPLETED
  - [x] `matrix_free::bicgstab`: Fix p_hat logic and test expectations
  - [x] `matrix_free::gmres`: Fix test expectations
  - [x] `multigrid::smoothers`: Fix Chebyshev eigenvalue estimation and bounds
  - [x] `spatial::weno`: Fix epsilon, test bounds, and test logic (negation removal, interface location)
  - [x] `time_stepping::imex`: Replace unstable ARK436L2SA with verified ARS343 (L-stable 3rd order)
  - [x] `time_stepping::stability`: Fix test assertion types
  - [x] `time_stepping::rk_chebyshev`: Updated with correct Verwer/Sommeijer recurrence logic; later Sprint 1.96.51 verification runs the focused RKC tests without ignored-case dependence

### Phase 3: Performance Benchmarking (Week 5-6)
- [ ] **Task 3.1**: Benchmarking infrastructure
  - [ ] Design benchmark configuration system
  - [ ] Implement timing and profiling utilities
  - [ ] Add memory usage tracking
- [ ] **Task 3.2**: Scaling analysis framework
  - [ ] Implement weak/strong scaling benchmarks
  - [ ] Add parallel efficiency metrics
  - [ ] Create scaling visualization tools
- [ ] **Task 3.3**: Production validation suite
  - [ ] Design production-scale test cases
  - [ ] Implement automated regression detection
  - [ ] Add performance alerting system

### Phase 4: Integration and Documentation (Week 7-8)
- [ ] **Task 4.1**: Framework integration
  - [ ] Integrate all components into cohesive validation suite
  - [ ] Add configuration management
  - [ ] Implement validation pipeline automation
- [ ] **Task 4.2**: Documentation and examples
  - [ ] Create comprehensive user documentation
  - [ ] Add tutorial examples for each feature
  - [ ] Generate API reference documentation
- [ ] **Task 4.3**: Final validation and testing
  - [ ] End-to-end validation of complete framework
  - [ ] Performance benchmarking of validation suite
  - [ ] Code review and quality assurance

## Quality Gates

### Code Quality
- [ ] **QG-001**: All code compiles without warnings (Passed for cfd-math)
- [ ] **QG-002**: Test coverage >85% for new validation code
- [ ] **QG-003**: Comprehensive documentation with examples
- [ ] **QG-004**: Code follows established patterns and idioms

### Validation Quality
- [ ] **QG-005**: MMS solutions verified against analytical results
- [ ] **QG-006**: Richardson extrapolation produces accurate convergence rates
- [ ] **QG-007**: Performance benchmarks show expected scaling behavior
- [ ] **QG-008**: Validation reports are clear and actionable

### Performance Requirements
- [ ] **QG-009**: Validation suite runs within reasonable time limits
- [ ] **QG-010**: Memory usage remains bounded for large problems
- [ ] **QG-011**: No performance regression in core CFD operations

## Risk Mitigation
- **Risk**: Complex geometry MMS implementation challenges
  - **Mitigation**: Start with simple geometries, build incrementally
- **Risk**: Richardson extrapolation numerical stability issues
  - **Mitigation**: Extensive testing with known analytical solutions
- **Risk**: Performance benchmarking overhead
  - **Mitigation**: Make benchmarking optional and configurable
- **Risk**: Accumulated Technical Debt
  - **Mitigation**: Emergency audit phase added to address critical gaps (Schwarz, ILU, Docs, Numerical)

## Sprint Burndown Tracking
- **Total Tasks**: 17
- **Completed**: 14
- **Remaining**: 3
- **Sprint Velocity**: 3.0 tasks/week (Phase 1 complete, Phase 2 complete, Phase 2b/c complete)

## Daily Standup Template
**Yesterday**: Updated RKC implementation with correct Verwer recurrence. Later Sprint 1.96.51 removed direct scalar conversion residue and verified focused RKC tests.
**Today**: Proceed to Phase 3 (Performance Benchmarking).
**Blockers**: None for the focused RKC scalar-provider slice; broader benchmarking remains pending.
**Next**: Task 3.1 Benchmarking infrastructure.

## Sprint Retrospective (End of Sprint)
**What went well?** Resolved critical stability issues in IMEX.
**What could be improved?** Broader time-stepping storage still needs the later Leto vector migration.
**Lessons learned?** Numerical schemes require exact coefficient verification.
**Action items for next sprint?** Continue provider migration for remaining time-stepping modules and vector storage.
