# Matrix-Free Operators and Krylov Solvers

For high-resolution flows, CFDrs prefers **matrix-free** evaluation: the
linear operator is a trait, and the solver allocates only the action
function `A·x`.  Matrix-free applications decouple memory footprint from
discretization size.

```rust
pub trait LinearOp<F: FloatElement> {
    fn apply(&self, x: &NdArray<F, Ix3>, y: &mut NdArray<F, Ix3>);
}
```

A Krylov solver asks only for `apply`, so it composes with any of the
three discretizations (Stencil, FEM, Spectral).

## Iterative Solvers

`cfd-math::krylov` exposes BiCGSTAB, GMRES, and IDR(s).  All three reach
the same convergence target at the same iteration count regardless of
which discretization backs them.

```rust
let krylov = BiCgStab::<f64>::new();
let solution = krylov.solve(&op, &rhs, &preconditioner, 1e-6, 5000)?;
```

## Examples Referenced by This Chapter

- [`matrix_free_demo`](examples/matrix_free_demo.md) — matrix-free
  operator application paired with a Krylov solver.

## Further Reading

- For spatial discretization surfaces, time integration, and adaptive
  time-stepping, see [Spectral, FEM, MUSCL, and Matrix-Free Methods](numerics_and_solvers.md).
- [`cfd-math` source](../../crates/cfd-math/src/)
