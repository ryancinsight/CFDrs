# Objectives
- Enforce exact Dirichlet handling and strict positivity of conductances
- Maintain quadratic-loss relinearization via `R_eff = R + 2k|Q|`
- Add residual tracking, refine SPD heuristics, and adaptively select CG/BiCGSTAB
- Expand manufactured-network and analytical validation tests
- Update rustdoc (invariants, boundary proofs, API) and CHANGELOG

# Implementation
## Boundary Conditions & Positivity
- Exact Dirichlet (row replacement + RHS substitution) is already implemented in `crates/cfd-1d/src/solver/matrix_assembly.rs:91–117`; retain and formalize invariants in module docs.
- Strict conductance positivity check exists in `crates/cfd-1d/src/solver/matrix_assembly.rs:72–77`; add explicit guard for `NaN/Inf` and epsilon floor to prevent division by zero in rare edge cases.
- Keep non-linear relinearization at the edge level: `r_eff = R + 2*k*|Q|` in `crates/cfd-1d/src/network/wrapper.rs:245–253`; additionally assert `r_eff > 0` before reciprocals.

## Quadratic Loss Relinearization
- Derivation: for `ΔP = R·Q + k·Q²`, linearize around current `Q_k` → effective resistance `R_eff = R + 2k|Q_k|`; already used in parallel edge iterator (`wrapper.rs:245–253`). Ensure consistent absolute value usage and monotonic positivity.

## Residual Tracking & Convergence
- Compute residual `||Ax − b||₂` each iteration (currently in `crates/cfd-1d/src/solver/mod.rs:325–338`); store per-iteration residuals in a `NetworkState` field and expose via API.
- Update `ConvergenceChecker` (`crates/cfd-1d/src/solver/convergence.rs`) to support dual criteria: solution change norm and residual norm, both below tolerance; add relative residual option.
- Provide structured tracing for residuals and iteration progress; keep zero-copy updates.

## SPD Heuristics & Adaptive Solver Selection
- Refine SPD detection (existing in `crates/cfd-1d/src/solver/mod.rs:282–311`) to:
  - Ignore identity Dirichlet rows in diagonal-dominance check
  - Enforce non-positive off-diagonals and positive diagonal for Laplacian rows
- Solver selection: use CG when SPD; otherwise BiCGSTAB (already in `crates/cfd-1d/src/solver/mod.rs:313–319`). Expose method choice in result metadata.
- Preconditioning: continue Diagonal Jacobi (`crates/cfd-1d/src/solver/linear_system.rs:93–121`), with division-by-zero protection.

# Tests
## Manufactured Networks
- Extend `crates/cfd-1d/tests/manufactured_network.rs`:
  - Series additivity (present: `test_series_additivity` at lines 8–50) → add variants with different Dirichlet placements
  - Parallel conductance addition (present: `test_parallel_additivity` at lines 52–101) → add asymmetric branch resistances and verify totals
  - Junction conservation (present: `test_conservation_at_junction` at lines 103–142) → add multi-outlet cases and numerical edge cases
  - Quadratic-loss edges: set `quad_coeff = k > 0`; verify convergence and correct `R_eff` application

## Analytical Validations
- Single pipe (Hagen–Poiseuille): verify `Q = ΔP / R` using `R` from `resistance/models/hagen_poiseuille.rs`
- Parallel/series analytical composites: compare computed flows vs analytical predictions
- Positivity failure tests: inject invalid conductance/resistance and assert `InvalidConfiguration`
- SPD heuristic tests: construct non-symmetric coupling to trigger BiCGSTAB; verify method selection

# Rustdoc & CHANGELOG
- Update `crates/cfd-1d/src/solver/mod.rs` Rustdoc invariants (`127–134`) to include explicit proofs that Dirichlet row replacement preserves SPD and discrete conservation.
- Add boundary-handling proof sketch in `crates/cfd-1d/src/solver/matrix_assembly.rs` (top module docs) referencing exact row/column handling.
- Ensure crate-level API references in `crates/cfd-1d/src/lib.rs:44–114` document new residual tracking and method-selection metadata.
- Append entries to root `CHANGELOG.md` describing: Dirichlet enforcement verification, positivity assertions, residual tracking exposure, SPD heuristic refinement, new manufactured/analytical tests, and docs updates.

# Acceptance Criteria
- All manufactured-network tests pass with tight tolerances (`≤ 1e-12` where applicable)
- Analytical validations match closed-form predictions within specified tolerances
- Residual histories are exposed in the result and available for tracing
- SPD heuristic selects CG for Laplacian-structured matrices and BiCGSTAB otherwise
- Rustdoc updated with formal statements, assumptions, and proofs; CHANGELOG includes precise changes and rationale

# Notes
- No behavioral masking: failures throw explicit errors, never clamped or hidden
- No placeholders: all invariants, proofs, and tests are complete, with documented assumptions and limitations