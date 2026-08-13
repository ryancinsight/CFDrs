# cfd-math

Numerical methods for the [CFDrs](https://github.com/ryancinsight/CFDrs) simulation
suite: the linear algebra, discretization, and time-integration layer that the 1D,
2D, and 3D solvers build on.

## What it provides

- **`linear_solver`, `iterative`** — CG, BiCGSTAB, GMRES, and preconditioners.
- **`multigrid`** — algebraic multigrid coarsening, cycles, and smoothers.
- **`sparse`** — sparse matrix storage and sparse matrix-vector products.
- **`time_stepping`** — explicit and implicit integrators, including
  Runge-Kutta-Chebyshev for mildly stiff problems.
- **`fd`, `fd_extensions`, `high_order`** — finite-difference stencils.
- **`quadrature_rules`, `interp`** — quadrature and interpolation.
- **`nonlinear_solver`, `optimization`, `statistics`, `diagnostics`** — supporting
  numerics.
- **`simd`** — architecture-conditional kernels dispatched at runtime.

All routines are generic over the scalar element traits from `eunomia`.

## Features

| Feature | Default | Effect |
|---|---|---|
| `gpu` | no | Enables `cfd-core/gpu` and GPU-backed metric collection |

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
