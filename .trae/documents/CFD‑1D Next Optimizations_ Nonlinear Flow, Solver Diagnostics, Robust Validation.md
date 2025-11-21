## Goals
- Enforce mathematically correct nonlinear losses in networks (ΔP = R·Q + k·Q²) via per‑iteration relinearization.
- Improve solver robustness and observability (residual tracking, SPD detection, adaptive method selection).
- Expand invariant tests and manufactured network cases (series/parallel laws, conservation under BCs, ill‑conditioning stress tests).
- Strengthen documentation of invariants and discrete BC proofs.

## Changes to Implement
### Nonlinear Flow Handling
- Add `quad_coeff` to edges (already added when applicable) and compute effective resistance each iteration: `R_eff = R + 2 k |Q|`.
- Ensure positivity constraints (`R > 0`, `k ≥ 0`), reject invalid parameters during assembly.

### Solver Diagnostics & Selection
- Add residual tracking (‖Ax − b‖ and solution delta norms) per iteration in `NetworkSolver::solve_network`.
- Detect SPD heuristically (symmetric structure + positive diagonal) and prefer CG; otherwise use BiCGSTAB.
- Keep CG as default, but add dynamic fallback if SPD test fails.

### Tests & Validation
- Add manufactured network tests:
  - Series additivity: `R_total = R1 + R2`.
  - Parallel additivity: `G_total = G1 + G2`.
  - Conservation at junctions with mixed Dirichlet/Neumann BCs.
  - Ill‑conditioned topologies (star/ladder) to verify convergence under positivity constraints.
- Expand rectangular duct correlation tests across aspect ratios; require `Re` explicitly.

### Documentation
- Update rustdoc in `matrix_assembly.rs` and `solver/mod.rs` with:
  - Invariants (units, positivity, applicability ranges).
  - Exact Dirichlet enforcement derivation and discrete Neumann consistency proof.
  - Nonlinear relinearization rationale and convergence criteria.

## Deliverables
- Correct nonlinear resistance handling integrated into assembly path.
- Solver loop with residual metrics, SPD detection, and robust method selection.
- New comprehensive tests passing under strict tolerances.
- Enhanced rustdoc documenting mathematical invariants and BC enforcement.

## Implementation Order
1) Add diagnostics (residual norms, SPD heuristic) and adaptive method selection.
2) Enforce positivity and parameter validation; keep relinearization in edge iterator.
3) Add manufactured tests and ill‑conditioning stress cases.
4) Update module rustdocs for invariants and BC proofs.
5) Run full test suite; iterate until all pass under strict tolerances.