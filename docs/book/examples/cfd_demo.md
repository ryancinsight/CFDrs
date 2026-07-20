# Example: cfd_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example cfd_demo`  
**Source**: [`examples/cfd_demo.rs`](../../examples/cfd_demo.rs)

## What This Example Demonstrates

Baseline sanity check covering four primitives of the CFDrs stack:

| Concept | API |
|---|---|
| 3D incompressible flow field | `FlowField::<f64>::new(nx, ny, nz)` |
| Divergence, vorticity, KE, enstrophy | `FlowOperations::{divergence, vorticity, kinetic_energy, enstrophy}` |
| Laminar/transitional/turbulent classification | `ReynoldsNumber::new(re, FlowGeometry::Pipe)` |
| Sparse matrix assembly | `SparseMatrixBuilder`, `IterativeLinearSolver` |

## Key Code Snippet

```rust
use cfd_core::physics::fluid_dynamics::{FlowField, FlowOperations};

let flow_field = FlowField::<f64>::new(32, 32, 32);
let div  = FlowOperations::divergence(&flow_field.velocity); // ≈ 0 for incompressible
let ke   = FlowOperations::kinetic_energy(&flow_field.velocity);

use cfd_core::physics::values::{FlowGeometry, ReynoldsNumber};
let re = ReynoldsNumber::new(2300.0, FlowGeometry::Pipe)?;
assert!(re.is_transitional()); // 2100 < Re < 4000 for pipe flow
```

## Physics Background

The **continuity equation** for an incompressible flow requires ∇·**u** = 0.
A healthy initial condition should give a divergence close to machine epsilon.
The **enstrophy** Ω = ∫½|**ω**|² dV is a global vorticity measure used for
turbulence diagnostics; for a quiescent initial field it should be zero.

## Book Chapter

[← CFDrs Architecture and Problem Setup](../foundations.md)

