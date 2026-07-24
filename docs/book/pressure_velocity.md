# Chapter 3 — Pressure-Velocity Coupling and Time Integration

CFDrs implements two families of pressure-velocity coupling algorithms:
**segregated** (SIMPLE/PIMPLE family) and **monolithic** (coupled solver).
Time integration covers both explicit and implicit strategies.

## Pressure-Velocity Coupling

### SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)

SIMPLE alternates momentum and pressure-correction sweeps until convergence:

```text
1. Solve momentum equation with current pressure: u* from A·u* = b − ∇p
2. Solve pressure-correction equation: ∇·(1/A · ∇p') = ∇·u*
3. Correct: p ← p + p', u ← u* − (1/A)·∇p'
4. Repeat until ‖∇·u‖ < ε
```

```rust
use cfd_core::solver::{SimpleSolver, SolverConfig};

let solver = SimpleSolver::new(SolverConfig { max_iter: 200, tolerance: 1e-5 });
let state = solver.solve(&mesh, &boundary_conditions, &initial_conditions)?;
```

### PIMPLE (SIMPLE with outer loops)

PIMPLE extends SIMPLE with outer corrector loops, enabling large Courant
numbers in transient simulations:

```rust
use cfd_core::solver::{PimpleSolver, PimpleConfig};

let pimple = PimpleSolver::new(PimpleConfig {
    n_outer_loops: 3,
    n_correctors: 2,
    ..Default::default()
});
```

### SIMPLEC (SIMPLE-Consistent)

SIMPLEC modifies the momentum under-relaxation to improve convergence
for mildly non-linear problems:

```rust
use cfd_core::solver::{SimplecSolver};
```

## Time Integration

### Explicit Schemes

For hyperbolic-dominated flows:

| Scheme | Order | CFL limit |
|--------|-------|-----------|
| Forward Euler | 1st | 1 |
| 4th-order Runge-Kutta | 4th | ~2.8 |
| Low-storage RK4 | 4th | ~2.8 |

```rust
use cfd_math::time_stepping::runge_kutta::RK4;

let rk4 = RK4::new(dt);
```

### Implicit Schemes

For diffusion-dominated and stiff problems:

| Scheme | Order | Stability |
|--------|-------|-----------|
| Crank-Nicolson | 2nd | A-stable |
| BDF2 | 2nd | A-stable |
| Runge-Kutta Chebyshev (RKC) | 2nd | Explicit, extended stability |

### Adaptive Time Stepping

`cfd-math` monitors the spectral radius and adjusts `dt` to maintain the
target CFL number:

```rust
use cfd_math::time_stepping::AdaptiveTimeStepper;

let ts = AdaptiveTimeStepper::new()
    .cfl_target(0.8)
    .dt_min(1e-8)
    .dt_max(0.01);
```

## Examples

- [simplec_pimple_demo](examples/simplec_pimple_demo.md) — SIMPLEC and PIMPLE coupling comparison
- [adaptive_time_stepping_demo](examples/adaptive_time_stepping_demo.md) — adaptive dt control
- [2d_heat_diffusion](examples/2d_heat_diffusion.md) — implicit heat equation
