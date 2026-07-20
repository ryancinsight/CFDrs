# Example: enhanced_cfd_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example enhanced_cfd_demo`  
**Source**: [`examples/enhanced_cfd_demo.rs`](../../examples/enhanced_cfd_demo.rs)

## What This Example Demonstrates

Advanced 2D incompressible Navier-Stokes around a circular cylinder at Re = 100 000
with LES turbulence, WENO5/7 reconstruction, and immersed-boundary geometry.

| Feature | API |
|---|---|
| Sub-grid scale models | `VremanModel`, `SigmaModel`, `MilesLES` |
| High-order reconstruction | `WenoReconstruction` (WENO5, WENO7) |
| Stiff time integration | `RungeKuttaChebyshev` |
| Immersed boundary | `ImmersedBoundaryMethod` |

## Key Code Snippet

```rust
use cfd_2d::physics::turbulence::{MilesLES, SigmaModel, VremanModel};
use cfd_math::high_order::weno::WenoReconstruction;
use cfd_math::time_stepping::RungeKuttaChebyshev;

let vreman = VremanModel::new();  // C_V coefficient from Vreman (2004)
let sigma  = SigmaModel::new();   // σ-model (Nicoud et al. 2011)
let miles  = MilesLES::new();     // MILES: monotone implicit LES

let weno = WenoReconstruction::weno5(); // 5th-order shock-capturing
let rk_cheb = RungeKuttaChebyshev::new(stages=10); // stable for stiff diffusion
```

## Physics Background

- **Vreman SGS**: `νt = C_V √(B_β/α_ij α_ij)`, zero in laminar regions by construction
- **WENO5**: 5th-order accurate for smooth regions, ENO near discontinuities
- **RK-Chebyshev**: s-stage explicit method stable up to `Δt ~ O(s²)` for diffusion

## Book Chapter

[← CFDrs Architecture and Problem Setup](../foundations.md)

