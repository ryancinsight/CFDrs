## Scope and Targets
- Crates: `cfd-core`, `cfd-math`, `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-mesh`, `cfd-io`, `cfd-validation`.
- Goals: Mathematical documentation, rigorous unit/property tests, physical-law validation, performance benchmarks, per-crate validation reports, and improved API docs.

## Phase 1: Audit & Model Inventory
- Enumerate operators, solvers, time-steppers, mesh/IO interfaces per crate.
- For each component, record: governing equations, discrete stencils, boundary conditions, invariants (symmetry/definiteness/nullspace), units, and assumptions.
- Identify gaps (mixed BC support, numerical consistency between CPU/GPU, missing invariants, absent tests).

## Phase 2: Test Design & Implementation
- Common test scaffolding (feature-gated GPU tests, proptest for properties, deterministic seeds).
- Floating-point precision analysis: absolute/relative error thresholds per operator; double-vs-single comparison where applicable; ulp-aware checks for critical paths.
- Dimensional analysis: runtime assertion checks (derived scalings), typed newtypes for spacing (dx, dy, dz) where feasible; consistent use of `1/dx²`, `1/dy²`, `1/dz²`.
- Conservation law verification: mass/momentum/energy where operators compose transport; divergence-free checks for incompressible flows; zero-mean constraints for Neumann/Periodic Laplacians.
- Numerical stability tests: time-stepping (explicit/implicit/IMEX) under CFL-like constraints; spectral radius/Von Neumann-style discrete amplification factors; regression tests for stability boundaries.
- Convergence tests: CG/GMRES/BiCGSTAB residual drops, theoretical bounds validation, SPD/Lanczos moment checks; multigrid cycle residual reduction factors and V/W-cycle behavior.
- Property-based tests (proptest): linearity, symmetry (<u, Lv>=<v, Lu>), scaling invariance, periodic wrap correctness, boundary stencil transitions; stochastic routine invariants (if any).

### Per-Crate Test Content
- `cfd-math`:
  - Operators: Laplacian (Dirichlet/Neumann/Periodic/mixed axes), Poisson3D, Momentum/Energy; MMS accuracy, invariants, GPU/CPU tolerance matching.
  - Solvers: CG/GMRES/BiCGSTAB — SPD tests, convergence-vs-tolerance, theoretical bounds, edge-case matrices.
  - Time stepping: stability regions, convergence order checks on manufactured ODEs/PDEs.
- `cfd-core`:
  - GPU binding correctness: uniform packing, bind group layout alignment with shaders, dispatch geometry; interface error types.
  - Buffer semantics and safety checks; basic CPU↔GPU round-trip validations (feature-gated).
- `cfd-1d/2d/3d`:
  - Domain-specific PDE stencils and BC suites; MMS across grids and spacings; conservation checks for transport.
  - Mixed boundary conditions (per-axis) and solver integration tests.
- `cfd-mesh`:
  - Geometry/topology: areas/volumes, adjacency, boundary indexing; partitioning consistency; refinement/coarsening invariants.
  - Dimensional consistency across coordinate transforms; I/O mapping to indices.
- `cfd-io`:
  - Schema validation, round-trip serialization fidelity (f32/f64), dimensional fields correctness, boundary/mesh metadata integrity.
- `cfd-validation`:
  - Aggregator harness to run crate suites, collect metrics, and emit structured validation summaries.

## Phase 3: Performance Benchmarks
- Criterion benches for computational kernels and solvers:
  - Operators: Laplacian/Poisson/Momentum/Energy CPU paths across grid sizes; optional GPU benches behind `gpu` feature.
  - Solvers: CG/GMRES/BiCGSTAB on canonical SPD/tridiagonal and Poisson-like systems with varying sizes.
- Run benches in `--release`; store results and trends per crate.

## Phase 4: Validation Reports
- Add `docs/VALIDATION_REPORT.md` per crate summarizing:
  - Mathematical models, discrete stencils, invariants, units
  - Test matrix with results (accuracy/stability/conservation/convergence)
  - Performance bench highlights
  - Assumptions, limitations, and known gaps
  - Cross-references to literature and APIs

## Phase 5: API Documentation
- Enhance rustdoc for public APIs with equations, examples, and invariants; add doc tests where safe.
- Ensure module-level docs for operator and solver traits, boundary condition semantics, and numerical properties.

## Phase 6: Tooling and CI Integration
- Prefer MSVC toolchain on Windows for test/linker reliability; gate GPU tests with `feature="gpu"`.
- Integrate `cargo nextest` and coverage (e.g., `tarpaulin` on Linux CI); enable per-crate test jobs and feature matrices.
- Bench jobs in CI scheduled or opt-in; upload artifacts for regression tracking.

## Deliverables
1. Per-crate `docs/VALIDATION_REPORT.md` with models, tests, benchmarks, and limitations.
2. Expanded test suites achieving full coverage of core algorithms (unit + property + integration + GPU-gated).
3. Criterion benchmarks and recorded performance results for kernels and solvers.
4. Enhanced API documentation with mathematical interfaces and examples.

## Execution Plan
- Iteration 1: Audit + implement tests for `cfd-math` and `cfd-core`; generate reports; run benches.
- Iteration 2: Extend tests/reports to `cfd-1d/2d/3d`; add conservation and MMS suites; time-stepping stability.
- Iteration 3: Mesh and IO validations; property tests for geometry/topology and serialization fidelity.
- Iteration 4: Validation aggregator in `cfd-validation` and CI/coverage integration.

## Notes
- Numerical thresholds will be documented per test; GPU/CPU tolerance windows defined by stencil order and FP precision.
- Mixed boundary axes standardized across operators; tests include periodic wrap and near-boundary stencil transitions.
- Where stochastic elements exist, property-based assertions will target invariants and distributional sanity checks.
