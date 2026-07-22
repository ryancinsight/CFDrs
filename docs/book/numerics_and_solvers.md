# Numerical Methods, Solvers, and Discretizations

CFDrs's solver layer decomposes into **discretization** (stencil/FEM),
**linear algebra** (matrix-free, sparse direct, Krylov), and **time
integration** (explicit, implicit, adaptive).  Commonalities across the
crates — including trait frontiers, error budgets, and atlas backends —
live in [`cfd-math`] so that they are not redeclared per solver.

## Spatial Discretization Surfaces

Three surfaces dominate CFDrs:

| Surface | Crate | Backend trait |
|---|---|---|
| Stencil (finite-volume) | `cfd-2d::stencil`, `cfd-3d::stencil` | `StencilOp<F>` |
| FEM (P1, P2 tetrahedral) | `cfd-3d::fem` | `FEMOp<F>` |
| Spectral (FFT-based) | `cfd-3d::spectral`, `cfd-2d::spectral` | `SpectralOp<F>` |

All three consume the boundary-and-state triple introduced in
[Foundations](foundations.md).

## Matrix-Free Operator Application

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
three discretizations.

## Iterative Solvers

`cfd-math::krylov` exposes BiCGSTAB, GMRES, and IDR(s).  All three reach
the same convergence target at the same iteration count regardless of
which discretization backs them.

```rust
let krylov = BiCgStab::<f64>::new();
let solution = krylov.solve(&op, &rhs, &preconditioner, 1e-6, 5000)?;
```

## Time Integration

```rust
pub trait Integrator<F: FloatElement> {
    fn step(&mut self, x: &mut NdArray<F, Ix3>, dt: F);
}

pub struct ForwardEuler<F>(PhantomData<F>);            // explicit
pub struct AdamsBashforth2<F>(PhantomData<F>);         // explicit, 2nd order
pub struct CrankNicolson<F>(PhantomData<F>);          // implicit
pub struct AdaptiveRungeKutta<F>(...);                 // adaptive
```

Adaptive Runge-Kutta is the default for stiff boundary-driven flows; the
explicit Adams-Bashforth pair is the default for incompressible time-stepping.

## Adaptive Time-Stepping

```rust
let stepping = AdaptiveStepping::builder()
    .cfl_max(0.5)
    .atol_field(1e-6)
    .build();

stepping.run(&mut solver, &state)?;
```

The adaptive scheme refines `dt` until local truncation error falls below
`atol_field` **per cell**.  This matches the Atlas Krylov trait frontier,
so callers can swap solvers without re-tuning tolerances.

## Examples Referenced by This Chapter

Part V opens with six example chapters:

- [`spectral_3d_poisson`](examples/spectral_3d_poisson.md) — spectral
  Poisson with FFT-plane decoupling.
- [`fem_3d_stokes`](examples/fem_3d_stokes.md) — 3-D Stokes FEM, P1 tetrahedral.
- [`matrix_free_demo`](examples/matrix_free_demo.md) — matrix-free operator
  application.
- [`muscl_schemes_demo`](examples/muscl_schemes_demo.md) — MUSCL limiters
  for the 2-D solver.
- [`simplec_pimple_demo`](examples/simplec_pimple_demo.md) — SIMPLEC and
  PIMPLE pressure-correction families.
- [`2d_heat_diffusion`](examples/2d_heat_diffusion.md) — explicit
  diffusion on a coarse grid.

## Further Reading

- [`cfd-math` source](../../crates/cfd-math/src/)
- [Geometry, Meshing, and CSG](geometry_and_meshing.md) for boundary
  generation upstream of these solvers.
