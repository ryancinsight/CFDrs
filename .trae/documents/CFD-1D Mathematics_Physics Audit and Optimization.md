## Summary of Current State
- 1D CFD is implemented as a steady-state network solver (no PDE discretization), assembling a sparse conductance matrix `A` and RHS `b` and solving `Ax=b`.
- Mathematical foundations are documented with circuit analogies and conservation laws; implementation aligns with linear laminar hydraulics across edges.

## Key Code References
- Boundary assembly: `crates/cfd-1d/src/solver/matrix_assembly.rs:59-69` (Dirichlet big-number penalty), `:70-75` (Neumann as source).
- Linear solver: `crates/cfd-1d/src/solver/linear_system.rs:60-85` (CG/BiCGSTAB, no preconditioner).
- Solver loop: `crates/cfd-1d/src/solver/mod.rs:259-300` (iterative update, convergence), invariants and SPD claims: `:156-164`.
- Flow update law: `crates/cfd-1d/src/network/wrapper.rs:220-229` (Q = (P_from − P_to)/R).
- Resistance models: Hagen–Poiseuille `crates/cfd-1d/src/resistance/models/hagen_poiseuille.rs:101-108`; Darcy–Weisbach `crates/cfd-1d/src/resistance/models/darcy_weisbach.rs:146-155` + friction factor `:172-238`; Rectangular duct `crates/cfd-1d/src/resistance/models/rectangular.rs:114-117` with Poiseuille number `:147-189`; Entrance effects `crates/cfd-1d/src/resistance/models/entrance.rs:125-133` (ΔP = R Q² mapping).

## Audit Findings
- Boundary enforcement uses a big-number diagonal penalty for Dirichlet (`1e20`). This is not mathematically sound across scales and can degrade conditioning, violating SPD assumptions.
- Neumann is handled as a source in `b`, which is consistent; however documentation does not prove the discrete consistency of source terms with conservation.
- Nonlinear loss modeling (entrance effects) defines resistance such that `ΔP = R Q²` but network flow update assumes `ΔP = R Q` linearly. Without relinearization of `R(Q)` per iteration, this is mathematically inconsistent.
- Rectangular channel uses Poiseuille correlations correctly; defaulting `Re = 100` when missing (`rectangular.rs:109-113`) hides applicability errors and can mask misuse.
- Solver claims SPD (`solver/mod.rs:162-164`), but SPD can be broken by mixed BC penalties or by negative/zero edge resistances; there are checks for zero resistance (`wrapper.rs:213-219`) but not for negative or scale-induced ill-conditioning.
- No preconditioning in iterative solvers, despite potential ill-conditioning from network topology and diameter scaling; this risks slow convergence or failure.
- Documentation is strong on physics background but lacks invariant statements for unit consistency, positivity, and discrete conservation proofs.
- Validation suite is solid for laminar hydraulics and correlations, but several tests are ignored (rectangular resistance behavior, series resistance, Colebrook convergence), indicating open correctness/performance gaps.

## Optimization & Correctness Objectives
1. Replace Dirichlet big-number method with mathematically sound enforcement:
   - Row-replacement (set row to identity and RHS to value) for exact Dirichlet.
   - Alternatively, Lagrange multipliers for constraints without conditioning damage.
2. Enforce rigorous invariants:
   - Units: ensure `R` units align with linear law throughout network code; prohibit quadratic `R(Q)` usage in edges unless treated as nonlinear.
   - Positivity: assert `R > 0`, `G = 1/R > 0` for all edges; forbid negatives.
   - Conservation: prove and test discrete mass conservation at all junctions for applied BCs.
3. Nonlinear components re-integration:
   - Introduce per-iteration relinearization: compute `R_eff = d(ΔP)/dQ` for components with `ΔP(Q)` and assemble with `R_eff` each iteration; update until convergence.
   - Provide a clear interface for components to supply `ΔP(Q)` and `dΔP/dQ`.
4. Solver robustness:
   - Add Jacobi/Block-Jacobi preconditioner for CG (SPD) and BiCGSTAB.
   - Auto-select CG when matrix is detected SPD; fall back to BiCGSTAB otherwise.
   - Add residual, conditioning metrics, and fail-fast diagnostics.
5. Boundary conditions rigor:
   - Prove and implement consistent Neumann handling (balanced RHS) and Robin if needed; document discrete derivations.
6. Rectangular duct correctness:
   - Require `Re` input; remove silent default; add applicability checks and errors when outside laminar range.
   - Fix ignored tests; cross-validate against Shah–London values across aspect ratios.
7. Validation expansion:
   - Network-level manufactured solutions: series/parallel networks with analytical solutions; ensure additivity laws.
   - Dimensional analysis tests: scaling with `L`, `D`, `ρ`, `μ` for all models.
   - Stress tests for ill-conditioned networks (star, ladder) to validate solver + preconditioning.
8. Documentation upgrades:
   - Add invariant section per module: Units, bounds, applicability ranges, discrete proofs.
   - Clarify differences between linear vs quadratic loss modeling and how they’re solved.
9. Performance:
   - CSR assembly micro-optimizations; avoid repeated pushes for BC rows by segregating BC application stage.
   - Optional parallel SPMV in solver config when beneficial.

## Implementation Plan (Phased)
- Phase A: Mathematical corrections
  1) Implement exact Dirichlet enforcement via row-replacement in `matrix_assembly.rs` and prove SPD preservation.
  2) Add strict invariant checks for edge resistances and unit consistency.
- Phase B: Nonlinear handling
  3) Extend edge/component API to support `ΔP(Q)` and `dΔP/dQ);` relinearize per iteration in `solve_network`.
- Phase C: Solver improvements
  4) Add Jacobi/Diagonal preconditioners and SPD detection; route to CG where applicable.
  5) Add robust convergence diagnostics and failure modes.
- Phase D: Model applicability & tests
  6) Require `Re` for rectangular model; remove fallback and enforce ranges.
  7) Un-ignore and fix failing tests; add manufactured-network tests for conservation/additivity.
- Phase E: Documentation
  8) Add module docs for invariants and discrete proofs; boundary condition derivations.
- Phase F: Performance
  9) Optimize assembly and optional parallel SPMV; benchmark network topologies.

## Deliverables
- Corrected BC enforcement and nonlinear component handling
- Preconditioned iterative solver with SPD-aware routing
- Complete invariant documentation and proofs
- Expanded test suite with no ignored tests and strict tolerances
- Benchmarks demonstrating improved convergence and stability