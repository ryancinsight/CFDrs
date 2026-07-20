# Example: turbulence_models_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example turbulence_models_demo`  
**Source**: [`examples/turbulence_models_demo.rs`](../../examples/turbulence_models_demo.rs)

## What This Example Demonstrates

Side-by-side demonstration of three RANS turbulence models on a 20×20 shear-flow grid,
showing how eddy viscosity and closure behaviour differ.

| Model | API | Regime |
|---|---|---|
| Standard k-ε (two-equation) | `KEpsilonModel` | Industrial internal flows |
| k-ω SST (Menter) | `KOmegaSSTModel` | Adverse-pressure-gradient flows |
| Spalart-Allmaras (one-equation) | `SpalartAllmaras` | Aerospace boundary layers |

## Key Code Snippet

```rust
use cfd_2d::physics::turbulence::{
    KEpsilonModel, KOmegaSSTModel, SpalartAllmaras, TurbulenceModel,
};
use cfd_2d::grid::StructuredGrid2D;

let grid = StructuredGrid2D::<f64>::new(20, 20, 0.0, 2.0, 0.0, 2.0)?;
// Each model is exercised over the same grid and initial conditions
// so eddy-viscosity magnitudes can be compared directly.
```

## Physics Background

All three models solve transport equations for turbulent scalars that feed into
an eddy-viscosity `νt` used in the momentum equation as `τ_ij = νt S_ij`.

- **k-ε**: robust for free shear; over-predicts separation
- **k-ω SST**: blends near-wall and free-stream; best for separated flows
- **SA**: minimal CPU cost; designed for attached aerodynamic flows

## Book Chapter

[← Turbulence Models and Cavitation](../turbulence_multiphase.md)

