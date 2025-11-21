## Objectives
- Reintroduce manufactured network tests using only public APIs to validate series/parallel additivity and junction conservation.
- Keep solver adaptive selection and residual diagnostics; avoid further intrusive changes.

## Implementation Steps
1) Add tests in `crates/cfd-1d/tests/manufactured_network.rs`:
   - Series network: two pipes in series; verify `Q = ΔP / (R1 + R2)` and equal branch flows.
   - Parallel network: two pipes from inlet to junction with low-resistance merge to outlet; verify `Q_total = ΔP (1/R1 + 1/R2)`.
   - Junction conservation: one inlet splitting to two outlets; verify inflow equals outflow at junction.
   - Use only public exports from `cfd-1d` (`NetworkBuilder`, `Edge`, `EdgeType`, `Network`, `EdgeProperties`, `NetworkProblem`, `NetworkSolver`, `network::ComponentType`) and `cfd_core::solver::Solver`.
2) Run the test suite and iterate until green.

## Deliverables
- New manufactured tests in `cfd-1d` validating additivity and conservation.
- No changes to public APIs or existing behavior; tests serve as mathematical verification.