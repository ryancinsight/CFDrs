# Example: bifurcation_2d_blood_validation

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example bifurcation_2d_blood_validation`  
**Source**: [`examples/bifurcation_2d_blood_validation.rs`](../../../examples/bifurcation_2d_blood_validation.rs)

## What This Example Demonstrates

Two 2D bifurcation validation cases (symmetric Casson, asymmetric Carreau-Yasuda)
solved on a staggered grid via the `BifurcationSolver2D` / SIMPLE path.

| Case | Rheology | Mass-balance limit | Symmetry limit |
|---|---|---|---|
| Symmetric | `CassonBlood` | 5 % | 5 % |
| Asymmetric | `CarreauYasudaBlood` | 10 % (cross-fidelity) | — |

## Key Code Snippet

```rust
use cfd_2d::solvers::bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
use cfd_2d::solvers::ns_fvm::{BloodModel, SIMPLEConfig};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};

const BLOOD_DENSITY: f64 = 1060.0;  // kg/m³
const INLET_VELOCITY: f64 = 0.1;   // m/s

let geom = BifurcationGeometry::symmetric_t_junction()?;
let casson = BloodModel::Casson(CassonBlood::default());
let result = BifurcationSolver2D::solve(&geom, &casson, INLET_VELOCITY, BLOOD_DENSITY,
                                        &SIMPLEConfig::default())?;

assert!(result.mass_balance_error < SYMMETRIC_MASS_BALANCE_LIMIT);
```

## Physics Background

At a symmetric T-junction, Murray's Law predicts equal flow in both daughters.
The **mass-balance error** ε = |Q_parent − Q_d1 − Q_d2| / Q_parent should be
< 5 % at validated resolution; the Cartesian staggered grid contributes ~5 %
geometric measurement error for angled branches.

## Book Chapter

[← Blood Flow and Rheology Workflows](../biomedical_flows.md)

