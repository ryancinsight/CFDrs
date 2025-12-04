# gap_audit

## Workflow
- identified → in_progress → resolved → validated (proof+test) → rejected

## Entries

### Unphysical Coefficients Accepted
- type: Error Masking
- status: validated (proof+test)
- context: Graph edges allowed `resistance=0` and `quad_coeff=0` simultaneously, or negatives.
- resolution: Early validation added in `NetworkBuilder::build` and `Network::validate_coefficients`; solver validates via `Validatable::validate_problem`.
- verification: Unit tests execute builder and solver on invalid configs expecting `InvalidConfiguration`. Release-mode tests pass.

### Quadratic Loss Linearization Consistency
- type: Math Errors
- status: validated (proof+test)
- context: Effective resistance for `ΔP = R·Q + k·Q|Q|` must be `R_eff = R + 2k|Q_k|`.
- resolution: Confirmed in `edges_parallel` and tested against closed-form per-branch flows.
- verification: `nonlinear_linearization_and_dirichlet.rs` validates `G = 1/R_eff`; manufactured parallel test validates branch flows via quadratic formula.

### Dirichlet Enforcement at Interior Junctions
- type: Algo Issues
- status: validated (proof+test)
- context: Identity row and RHS injection must hold for interior fixed nodes; neighbors must accumulate diagonal `ΣG` and RHS `ΣG·P_dir`.
- resolution: Assembly enforces exact Dirichlet; new test covers interior junction.
- verification: Tests assert identity row and correct RHS/diagonal for neighbors.

### Solver Method Selection Correctness
- type: Working-But-Incorrect
- status: resolved
- context: Ensure CG only on SPD active submatrix; otherwise BiCGSTAB.
- resolution: Heuristic SPD check excludes Dirichlet identity rows and enforces Laplacian structure.
- verification: Existing tests and residual tracking confirm stable selection; no failures.

