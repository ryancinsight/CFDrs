# Example: SIMPLEC and PIMPLE Pressure-Velocity Coupling

**Crate**: `cfd-suite` (workspace root)
**Run**: `cargo run --example simplec_pimple_demo`
**Source**: [`examples/simplec_pimple_demo.rs`](../../../examples/simplec_pimple_demo.rs)

## What This Example Demonstrates

Demonstrates SIMPLEC (consistent SIMPLE) and PIMPLE (merged PISO-SIMPLE) for incompressible flow.

| API | Purpose |
|---|---|
| `cfd_2d::solvers::{SimplecSolver, PimpleSolver}` | SIMPLEC (Van Doormaal & Raithby 1984): consistent Rhie-Chow interpolation |
| `SIMPLEConfig` | SIMPLEC (Van Doormaal & Raithby 1984): consistent Rhie-Chow interpolation |

## Physics Background

SIMPLEC (Van Doormaal & Raithby 1984): consistent Rhie-Chow interpolation. PIMPLE: outer PISO correctors + inner SIMPLE for transient stiff flows.

## Book Chapter

[← Part V — Discretization and Solvers](../numerics_and_solvers.md)