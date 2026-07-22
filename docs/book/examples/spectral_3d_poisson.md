# Example: spectral_3d_poisson

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example spectral_3d_poisson`  
**Source**: [`examples/spectral_3d_poisson.rs`](../../../examples/spectral_3d_poisson.rs)

## What This Example Demonstrates

Solves ∇²u = f on a 3D unit cube with Dirichlet boundary conditions using the
spectral Poisson solver, and compares against the analytical solution.

| Concept | API |
|---|---|
| 3D Poisson solver creation | `PoissonSolver::new(nx, ny, nz)` |
| Dirichlet boundary specification | `PoissonBoundaryCondition::Dirichlet` |
| RHS and solution arrays | `leto::Array1::zeros([n])` |

## Key Code Snippet

```rust
use cfd_3d::spectral::{PoissonBoundaryCondition, PoissonSolver};
use leto::Array1;
use std::f64::consts::PI;

let solver = PoissonSolver::new(8, 8, 8)?;

// f(x,y,z) = sin(πx)sin(πy)sin(πz)
// Exact: u = -sin(πx)sin(πy)sin(πz) / (3π²)
let mut rhs = Array1::zeros([8 * 8 * 8]);
for i in 0..8 { for j in 0..8 { for k in 0..8 {
    let x = -1.0 + 2.0 * i as f64 / 7.0;
    let y = -1.0 + 2.0 * j as f64 / 7.0;
    let z = -1.0 + 2.0 * k as f64 / 7.0;
    rhs[i*64 + j*8 + k] = (PI*x).sin() * (PI*y).sin() * (PI*z).sin();
}}}
```

## Physics Background

The **spectral Poisson solver** represents u as a truncated Fourier series and
solves in wavenumber space via division by −|**k**|². This achieves spectral
convergence (exponential error decay) for smooth f, which makes it the
preferred baseline for pressure-Poisson steps in incompressible CFD.

Zero padding: the solver uses Atlas `apollo-fft` for the forward/inverse
transforms so the convolution is computed exactly without aliasing.

## Book Chapter

[← Spectral, FEM, MUSCL, and Matrix-Free Methods](../numerics_and_solvers.md)

