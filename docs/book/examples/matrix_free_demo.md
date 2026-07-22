# Example: matrix_free_demo

**Crate**: `cfd-suite` (workspace root)  
**Run**: `cargo run --example matrix_free_demo`  
**Source**: [`examples/matrix_free_demo.rs`](../../../examples/matrix_free_demo.rs)

## What This Example Demonstrates

Matrix-free CG and Laplacian operator application using the `LinearOperator`
trait, solving a 1D diffusion problem and a 2D Laplacian system without ever
storing the coefficient matrix.

| API | Purpose |
|---|---|
| `LinearOperator<f64>` trait | Implement a matrix-free operator |
| `LaplacianOperator2D` | Typed provider-backed negative 2D Laplacian |
| `ConjugateGradient` | Iterative solver for SPD systems |
| `IterativeSolverConfig` | Tolerance, max-iteration control |

## Key Code Snippet

```rust
use cfd_math::linear_solver::{
    ConjugateGradient, IterativeSolverConfig, LaplacianOperator2D, LinearOperator,
};
use leto::Array1;

// 1D diffusion operator (d²u/dx²) with homogeneous Dirichlet BCs
struct DiffusionOperator1D { n: usize, dx: f64 }

impl LinearOperator<f64> for DiffusionOperator1D {
    fn apply(&self, x: &Array1<f64>, y: &mut Array1<f64>) -> Result<()> {
        let dx2_inv = 1.0 / (self.dx * self.dx);
        for i in 1..(self.n - 1) {
            y[i] = (x[i-1] + x[i+1] - 2.0*x[i]) * dx2_inv;
        }
        y[0] = 0.0; y[self.n-1] = 0.0;
        Ok(())
    }
    fn size(&self) -> usize { self.n }
}

let cg = ConjugateGradient::new(IterativeSolverConfig::default());
let sol = cg.solve(&op, &rhs)?;
```

## Why Matrix-Free?

For large structured grids the stencil can be applied in O(N) memory and time
without storing an N×N matrix. The `LinearOperator` trait enables any
operator — spectral, FEM, or finite-difference — to plug into the same CG/GMRES
iterative solver path.

## Book Chapter

[← Spectral, FEM, MUSCL, and Matrix-Free Methods](../numerics_and_solvers.md)

