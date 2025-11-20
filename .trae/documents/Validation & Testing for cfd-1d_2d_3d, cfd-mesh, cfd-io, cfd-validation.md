## Targets and Goals
- Crates: `cfd-1d`, `cfd-2d`, `cfd-3d`, `cfd-mesh`, `cfd-io`, `cfd-validation`.
- Deliverables: Mathematical documentation (validation reports), rigorous unit/property/integration tests, physical-law validation, performance benchmarks (where applicable), and an aggregated validation harness.

## Per-Crate Work
### cfd-1d / cfd-2d / cfd-3d
- Document PDEs, stencils, boundary conditions, invariants, and units in `docs/VALIDATION_REPORT.md`.
- Implement MMS accuracy tests for representative operators (Laplacian/Poisson/momentum/energy) with Dirichlet/Neumann/Periodic/mixed boundaries.
- Convergence tests for iterative solvers integrated with these operators; verify theoretical residual bounds.
- Conservation checks for transport (mass/momentum/energy) where applicable.
- GPU feature-gated tests matching CPU math with tolerance thresholds (if GPU paths exist).
- Criterion benches for operator apply and small solver runs in `--release`.

### cfd-mesh
- Tests for geometry/topology correctness: cell areas/volumes, boundary indexing, adjacency consistency.
- Dimensional analysis and unit checks across coordinate transforms and mesh refinement/coarsening operations.
- Validation report summarizing models, invariants, tests, and assumptions.

### cfd-io
- Round-trip serialization/deserialization tests (f32/f64) ensuring no precision loss beyond documented thresholds.
- Schema validation tests for boundary/mesh metadata; dimensional field checks.
- Validation report documenting assumptions, formats, and invariants.

### cfd-validation
- Create a validation harness to orchestrate running per-crate test suites and collect summaries.
- Emit consolidated validation report (table/summary with coverage and key metrics), reusing per-crate reports.

## Testing Requirements
- Floating-point precision analysis with explicit absolute/relative tolerances per operator and boundary mode.
- Dimensional analysis checks for all physical computations.
- Conservation law verification for transport equations.
- Numerical stability tests for time-stepping interfaces; CFL-style constraints and robustness checks.
- Property-based testing (linearity, symmetry, scaling, periodic wrap, boundary transitions; stochastic invariants where present).

## Performance Benchmarks
- Criterion benchmarks for operators and solvers across grid sizes; record timings and trends.

## Execution Order
1. Implement `cfd-1d/2d/3d` validation reports, MMS tests, and solver integrations.
2. Add `cfd-mesh` geometry/topology tests and report.
3. Add `cfd-io` serialization and schema tests and report.
4. Implement `cfd-validation` aggregator and consolidated report.

## Notes
- Use feature gating for GPU tests; maintain MSVC toolchain guidance on Windows for reliable linking.
- Document tolerances and limitations explicitly in reports; ensure API docs reflect mathematical interfaces and invariants.