# Chapter 2 — Governing Equations and Discretization

CFDrs solves the incompressible Navier-Stokes equations in primitive-variable
form together with optional scalar transport equations (energy, turbulent
kinetic energy, species concentration).

## Incompressible Navier-Stokes

The continuity and momentum equations for a constant-density fluid:

```text
∇·u = 0

∂u/∂t + (u·∇)u = -∇(p/ρ) + ν ∇²u + f/ρ
```

where **u** is velocity, *p* pressure, *ρ* density, *ν* kinematic viscosity,
and **f** body forces.  CFDrs resolves the incompressibility constraint through
the **pressure-Poisson equation** (see [Chapter 3](pressure_velocity.md)).

## Energy Equation

For thermal flows, temperature *T* is transported:

```text
∂T/∂t + (u·∇)T = α ∇²T + S_T
```

with thermal diffusivity α = λ/(ρ c_p).  Density-driven buoyancy is
handled through the Boussinesq approximation: f_y = -ρ g β (T - T_ref).

## Non-Newtonian Rheology

Blood and polymer flows use the generalised Newtonian model:

```text
τ = μ_eff(γ̇) · γ̇
```

`cfd-core` provides the rheology dispatch:

```rust
use cfd_core::physics::rheology::{Rheology, PowerLaw, Casson, CarreauYasuda};

let blood = Rheology::Casson { yield_stress: 0.015, mu_inf: 0.0035 };
```

Viscosity depends on shear rate γ̇ = |2**D**|, where **D** is the symmetric
velocity gradient tensor.

## Spatial Discretisation

CFDrs supports three spatial discretisation families, each owning a dedicated
crate:

| Method | Crate | Order |
|--------|-------|-------|
| Finite Volume (staggered / collocated) | `cfd-2d`, `cfd-3d` | 2nd |
| Spectral (Fourier / Chebyshev) | `cfd-3d::spectral` | spectral |
| FEM (P1, P2 tetrahedral) | `cfd-3d::fem` | 1st–2nd |
| MUSCL reconstruction | `cfd-3d::muscl` | up to 5th |

Boundary conditions live in `cfd-core::physics::boundary` and are consumed
by every discretisation layer through a single `BoundaryCondition<T>` enum:

```rust
use cfd_core::physics::boundary::{BoundaryCondition, BoundaryConditionSet};

let bc = BoundaryConditionSet::<f64>::new()
    .add("inlet",  BoundaryCondition::velocity_inlet([1.0, 0.0, 0.0].into()))
    .add("walls",  BoundaryCondition::wall_no_slip())
    .add("outlet", BoundaryCondition::pressure_outlet(0.0));
```

## Examples

- [cfd_demo](examples/cfd_demo.md) — end-to-end setup from grid to result
- [enhanced_cfd_demo](examples/enhanced_cfd_demo.md) — multi-physics setup with plugins
