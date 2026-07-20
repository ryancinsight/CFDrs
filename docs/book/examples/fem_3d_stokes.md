# Example: fem_3d_stokes

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example fem_3d_stokes`  
**Source**: [`examples/fem_3d_stokes.rs`](../../examples/fem_3d_stokes.rs)

## What This Example Demonstrates

Solves the creeping (Stokes) flow equations in a unit-cube lid-driven cavity
using P1/P1 stabilised tetrahedral finite elements.

| Concept | API |
|---|---|
| Constant-property fluid | `ConstantPropertyFluid::<f64>::water_20c()` |
| Cube mesh via Kuhn triangulation | `create_refined_cube_mesh(nx)` |
| FEM problem + solver | `StokesFlowProblem`, `FemSolver`, `FemConfig` |
| Boundary conditions | `BoundaryCondition::NoSlip`, `BoundaryCondition::Dirichlet` |

## Key Code Snippet

```rust
use cfd_3d::fem::{FemConfig, FemSolver, StokesFlowProblem};
use cfd_core::physics::fluid::ConstantPropertyFluid;

let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
// ρ = 998.2 kg/m³,  μ = 1.002e-3 Pa·s
println!("ν = {:.2e} m²/s", fluid.viscosity / fluid.density);

let mesh = create_refined_cube_mesh(5)?; // 5×5×5 → 125 verts, ~500 tets
let problem = StokesFlowProblem::new(&mesh, &fluid)?;
let solution = FemSolver::solve(&problem, &FemConfig::default())?;
```

## Physics Background

The **Stokes equations** are the zero-inertia limit of Navier-Stokes:

```
−∇p + μ ∇²u = 0,   ∇·u = 0
```

These are linear, making FEM solution direct (one linear system). They are the
correct model when Re ≪ 1 (microfluidics, polymer processing, geophysical
creeping flows). Reference geometry: Ghia et al. (1982).

## Book Chapter

[← Spectral, FEM, MUSCL, and Matrix-Free Methods](../numerics_and_solvers.md)

