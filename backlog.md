# CFDrs Backlog

## CURRENT PHASE: CFD-1D AUDIT (Sprint 1.100.0)
- **Goal**: Audit, optimize performance, enhance memory efficiency, correct, complete and attempt to prove implementation wrong in `cfd-1d`.
- **Status**: ✅ COMPLETE
- **Methodology**: 
  - Trace code paths for mathematical correctness under the Navier-Stokes / Hagen-Poiseuille 1D assumptions.
  - Review network solvers for memory bottlenecks (unnecessary `clone()`, `HashMap` vs `Vec` indexing).
  - Ensure strict physical constraint invariants are checked (e.g. Reynolds number limits, compressibility limits).
  - Write negative tests.

### Strategy Outline:
1. Validated and converted Linear System Solver efficiency (O(1) Boundary vector mappings).
2. Validated First-principles physical equations bounds.
3. Implemented Adverse/Boundary testing preventing silent `NaN` flows.
