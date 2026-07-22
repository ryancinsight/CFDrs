# Chapter 1 — CFDrs Architecture and Problem Setup

## Overview

CFDrs is a workspace of CFD-focused crates organized so that *physics,
geometry, and solver* concerns live in **distinct** crates and compose through
workspace-level contracts. Unlike monolithic CFD frameworks that entangle
mesh generation, discretization, linear algebra, and post-processing, CFDrs
adopts a strict layered architecture where each crate owns one bounded context
and re-exports no foreign implementation detail.

The central invariant is the **canonical problem triple**: `grid + boundaries +
state`. Every solver entry point consumes exactly these three inputs,
produces a converged state, and hands it to a posterior that computes
quantities of interest. This chapter formalizes that triple, introduces the
governing equations, taxonomizes boundary conditions, presents the simulation
lifecycle with its Atlas stack integration points, and records the
nondimensionalization conventions that govern every validation oracle.

---

## Governing Equations — Incompressible Navier-Stokes

CFDrs is centered on the incompressible Navier-Stokes equations in primitive
variables. The strong form on a domain $\Omega \subset \mathbb{R}^d$ with
boundary $\partial\Omega$ and time interval $[0, T]$ reads

$$\nabla \cdot \mathbf{u} = 0 \tag{1}$$

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f} \tag{2}$$

where $\mathbf{u}(\mathbf{x}, t) \in \mathbb{R}^d$ is velocity,
$p(\mathbf{x}, t) \in \mathbb{R}$ is pressure, $\rho$ is density (constant for
incompressible flow), $\nu = \mu / \rho$ is kinematic viscosity, and
$\mathbf{f}$ is a volumetric body force. Equation (1) is the incompressibility
constraint (divergence-free condition); equation (2) is the momentum balance
including convective transport $(\mathbf{u} \cdot \nabla)\mathbf{u}$,
pressure gradient, viscous diffusion $\nu\nabla^2\mathbf{u}$, and forcing.

### Vorticity-Streamfunction Form (2-D)

For 2-D incompressible flows CFDrs often uses the vorticity-streamfunction
formulation, which eliminates pressure as an explicit unknown. Define vorticity
$\omega = \nabla \times \mathbf{u} = \partial_x v - \partial_y u$ and
streamfunction $\psi$ satisfying $\nabla^2 \psi = -\omega$ with
$\mathbf{u} = (\partial_y \psi, -\partial_x \psi)$. The coupled system becomes

$$\frac{\partial \omega}{\partial t} + (\mathbf{u} \cdot \nabla)\omega = \nu \nabla^2 \omega \tag{3}$$

$$\nabla^2 \psi = -\omega \tag{4}$$

Equation (3) is a scalar advection-diffusion equation for $\omega$; equation (4)
is a Poisson equation for $\psi$. This formulation is exploited by
`cfd-2d::physics::vorticity_stream` and is the workhorse for lid-driven cavity
validation against Ghia et al. (1982).

### Reynolds Number

The Reynolds number

$$\mathrm{Re} = \frac{UL}{\nu} = \frac{\rho U L}{\mu} \tag{5}$$

with characteristic velocity $U$ and length $L$ parametrizes the ratio of
inertial to viscous forces. CFDrs canonical benchmarks sweep
$\mathrm{Re} \in [100, 10000]$ for cavity flows, $[100, 2100]$ for pipe flows
(laminar regime), and $[180, 590]$ (friction Reynolds $Re_\tau$) for turbulent
channel references.

---

## Nondimensionalization

All internal solvers operate in nondimensional variables unless a crate-level
doc comment states otherwise. The nondimensional scheme is:

$$\mathbf{x}^* = \frac{\mathbf{x}}{L}, \quad t^* = \frac{t\, U}{L}, \quad \mathbf{u}^* = \frac{\mathbf{u}}{U}, \quad p^* = \frac{p}{\rho U^2} \tag{6}$$

Substituting into (1)-(2) and dropping the $*$ yields

$$\nabla \cdot \mathbf{u} = 0, \qquad \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \frac{1}{\mathrm{Re}}\nabla^2\mathbf{u} + \mathbf{f}^* \tag{7}$$

For biomedical workflows where physical units matter (stenosis screening,
hemolysis), `cfd-1d` and `cfd-schematics` carry an explicit length scale and
convert to nondimensional variables only at the solver interface boundary via

```rust
pub struct NonDimParams<F> {
    pub re: F,       // Reynolds number
    pub length: F,   // characteristic length [m]
    pub velocity: F, // characteristic velocity [m/s]
}
```

Validation oracles always document which scaling they use.

---

## Boundary Condition Taxonomy

CFDrs models boundary conditions as a closed `enum` in `cfd-core::boundary`.
Every face or surface element in the mesh carries exactly one
`BoundaryKind<F>` value. The taxonomy is:

### Dirichlet (Prescribed Value)

- **`Wall`** — no-slip, no-penetration: $\mathbf{u} = \mathbf{0}$ on
  $\partial\Omega_{\mathrm{wall}}$.
- **`Lid { velocity }`** — moving wall (lid-driven cavity): $\mathbf{u} =
  (U_{\mathrm{lid}}, 0)$ with $U_{\mathrm{lid}}$ typically 1.0 in nondim
  units.
- **`Inlet { velocity_profile }`** — prescribed inflow, uniform or parabolic:
  $\mathbf{u} = \mathbf{u}_{\mathrm{in}}(y,z)$. For Hagen-Poiseuille pipe, the
  parabolic profile $u(r) = 2 U_{\mathrm{avg}}(1 - r^2/R^2)$ is imposed.
- **`InletPressure { p }`** — pressure-driven inflow.

In vorticity-streamfunction form the Dirichlet velocity BC translates to
$\omega_{\mathrm{wall}} = -2\psi_{i,wall}/\Delta n^2$ (Thom's formula, first
order) or higher-order Briley form.

### Neumann (Prescribed Flux/Gradient)

- **`Outflow`** — zero-gradient outflow: $\partial\mathbf{u}/\partial n = 0$.
- **`OutflowPressure { p }`** — fixed pressure at outlet: $p = p_{\mathrm{out}}$.
- **`Symmetry`** — symmetry plane: $\mathbf{u}\cdot\mathbf{n} = 0$,
  $\partial(\mathbf{u}\cdot\mathbf{t})/\partial n = 0$.
- **`TractionFree`** — zero traction: $\boldsymbol{\sigma}\cdot\mathbf{n} = 0$
  for FEM assembly.

### Robin (Mixed)

- **`PartialSlip { slip_length }`** — Navier slip: $\mathbf{u}\cdot\mathbf{t} =
  L_s \, \partial(\mathbf{u}\cdot\mathbf{t})/\partial n$ where $L_s$ is the slip
  length, relevant for hydrophobic microfluidics and porous scaffolds.
- **`PressureResistance { resistance }`** — lumped-parameter outlet:
  $p = R \, Q$ linking pressure to flowrate through a resistance.

### Time-Dependent

- **`Pulsatile { waveform }`** — Womersley inlet for biomedical flows:
  $\mathbf{u}(t) = \mathbf{u}_{\mathrm{mean}} + \sum_k
  \mathbf{A}_k e^{i(k\omega_0 t)}$. Used by `cfd-1d::WomersleyInlet`.

Each BC variant carries typed parameters (no bare `f64` domain values) and
implements `BoundaryCondition<F: FloatElement>` with methods `apply_dirichlet`
/ `apply_neumann` that the solver calls during matrix assembly.

---

## Crate Map

| Crate | Role | Key types |
|---|---|---|
| `cfd-core` | Scalar seams, dimensional types, boundary conditions, error types | `FloatElement`, `BoundaryKind`, `CfdError`, `StateVec` |
| `cfd-math` | Linear algebra, sparse matrices, iterative solvers, FFT bridge | `CsrMatrix`, `NdArray`, `BiCgStab`, `FftPlan` |
| `cfd-1d` | 1-D blood/Reynolds models, microfluidic screens, hemolysis pipelines | `NetworkSolver`, `CarreauModel`, `HemolysisIndex` |
| `cfd-2d` | 2-D Stencil solvers, biomass CFD, CSG-coupled channels | `StructuredGrid2D`, `VorticityStreamSolver`, `MusclLimiter` |
| `cfd-3d` | FEM and spectral 3-D solvers, multi-region flow | `TetMesh`, `SpectralOp`, `TurbulenceModel` |
| `cfd-validation` | Reference phantoms and orthogonal oracles | `GhiaOracle`, `PoiseuilleOracle`, `RichardsonExtrapolation` |
| `cfd-schematics` | CSG primitives and curve generators | `NetworkBlueprint`, `VenturiSpec`, `BifurcationBuilder` |
| `cfd-io` | HDF5 / NetCDF / CSV output adapters | `VtkWriter`, `CsvSink` |
| `cfd-optim` | Optimization audit harness | `DesignSpace`, `ParetoFront` |
| `cfd-python` | PyO3 bindings | thin `#[pyclass]` wrappers |

Every solver crate depends on `cfd-core` and `cfd-math`; geometry crates
(`cfd-schematics`) are depended on by `cfd-{1,2,3}d`. The cross-cutting
contracts — boundary conditions, time-stepping, state vectors — live in
`cfd-core` and are referenced by example chapters throughout this book.

---

## State Vector and Boundary Types

A CFDrs solver entry point looks like:

```rust
use cfd_core::boundary::{Boundary, BoundaryKind};
use cfd_core::state::StateVec;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::vorticity_stream::{VorticityStreamConfig, VorticityStreamSolver};

let grid = StructuredGrid2D::<f64>::unit_square(41, 41)?;

let mut state = StateVec::<f64, Ix2>::with_capacity(&grid)?;
state.set_pressure(1.0)?;

let boundaries = [
    Boundary::new(BoundaryKind::Lid { velocity: 1.0 }),
    Boundary::new(BoundaryKind::Wall),
    Boundary::new(BoundaryKind::Wall),
    Boundary::new(BoundaryKind::Wall),
];

let mut solver = VorticityStreamSolver::new(grid, boundaries, state)?;
solver.run_until_steady(1e-6, 10000)?;
```

The `grid + boundaries + state` triple is the **canonical problem setup**.
Every part of this book references one of those three components. `StateVec`
internally holds `leto::NdArray` (Atlas) or `ndarray::Array` (legacy) fields
for velocity components, pressure, vorticity, and optional auxiliary fields
(temperature, species mass fractions, phase fraction).

Validation workflow:

```rust
// Validate the converged field against Ghia et al.
let ghia = GhiaOracle::re_400();
let err = ghia.l2_error(&solver.velocity_profile(Centerline::Horizontal));
assert!(err < 1e-3, "Ghia deviation {err}");
```

---

## The CFDrs Simulation Lifecycle

```
+------------+   +-----------+   +----------+   +----------+
| Geometry   |-->| Mesh/Grid |-->| Boundary |-->| Initial  |
| (CSG)      |   |           |   | wiring   |   | state    |
+------------+   +-----------+   +----------+   +----------+
                                                    |
                                                    v
+------------+   +-----------+   +----------+   +----------+
| Posterior  |<--| Solve     |<--| Time     |   | Adapt    |
| analysis   |   | (linear   |   | step     |-->| (mesh /  |
|            |   |  iter.)   |   |          |   |  dt)     |
+------------+   +-----------+   +----------+   +----------+
       |                    ^                        |
       +--------------------+                        |
       |  validate against oracles                   |
       +---------------------------------------------+
```

Concretely:

1. **Geometry → Mesh** is handled in `cfd-schematics` (CSG) and
   `cfd-{2,3}d/grid.rs`. Primitives are combined via union/intersection/
   difference, then tessellated into a `Mesh` or `StructuredGrid2D/3D`.
2. **Boundary wiring** lives in [`cfd_core::boundary`]. Each mesh face gets a
   `BoundaryKind`. Higher-level helpers (`BifurcationBuilder::wire_boundaries`)
   automate wiring for common topologies.
3. **Initial state** is a `StateVec` filled with uniform, parabolic, or
   rest-state data.
4. **Time stepping** is configurable per solver (`cfd-2d::time` /
   `cfd-3d::time` / adaptive in `cfd-2d::adaptive`). Explicit Adams-Bashforth
   for advection, implicit Crank-Nicolson for diffusion, SIMPLEC/PIMPLE for
   pressure-velocity coupling.
5. **Solve** per step means Poisson solve for pressure or streamfunction (via
   `cfd-math::krylov`), transport solve for vorticity/velocity, turbulence
   closure update.
6. **Adapt** — mesh refinement (gradient-based or residual-based) and
   CFL-driven timestep adaptation. Richardson extrapolation quantifies spatial
   convergence order.
7. **Posterior** means DVH-equivalent (cavitation risk, shear envelope,
   Reynolds stress, hemolysis index), reported through user-supplied hooks.
8. **Validate** — compare converged fields against analytical or published
   numerical oracles (`cfd-validation`).

### Atlas Stack Interface

```
User code
  │
  ├─ cfd-core   ──► eunomia (traits) + mnemosyne (arena)
  ├─ cfd-math   ──► leto (NdArray, Csr) + hermes-simd + apollo (FFT)
  ├─ cfd-1d/2d/3d ► leto + hermes-simd + moirai (parallel) + apollo
  ├─ cfd-schematics ► leto::Point3/Vector3 + hephaestus (GPU opt-in)
  └─ cfd-io     ──► consus (planned)
```

Where a crate already declares Atlas crates in `Cargo.toml`, the migration is
**complete** for that surface; where `ndarray` or `nalgebra` remains in the
manifest, the crate is in **bulk-migration** phase. See [Atlas Dependency
Map](appendix_dependencies.md) for the per-crate table and
[Atlas Stack Integration — Migration Reference](performance_and_atlas.md) for
the migration chapters. Migration progress also lives in
[`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md) and the per-crate roadmap sections.

---

## Design Contract Enumeration

| Contract | Crate | Purpose |
|---|---|---|
| `FloatElement` | `cfd-core` / `eunomia` | Generic scalar bound (`f32`/`f64`) |
| `BoundaryKind<F>` | `cfd-core` | BC taxonomy |
| `StateVec<F,D>` | `cfd-core` | Field storage |
| `LinearOp<F>` | `cfd-math` | Matrix-free operator action |
| `TurbulenceModel` | `cfd-3d` | Turbulent viscosity closure |
| `Cavitation` | `cfd-3d` | Cavitation inception/collapse |
| `MeshGenerator<F>` | `cfd-schematics` | Mesh construction |
| `GeometryPipeline<F>` | `cfd-schematics` | Full CSG→solver pipeline |
| `Rheology` | `cfd-1d` | Non-Newtonian viscosity |

---

## Examples Referenced by This Chapter

- [`cfd_demo`](examples/cfd_demo.md) — minimum end-to-end setup.
- [`enhanced_cfd_demo`](examples/enhanced_cfd_demo.md) — diagnostics enabled.
- [`adaptive_time_stepping_demo`](examples/adaptive_time_stepping_demo.md)
  — demonstrates the time-stepping hook with adaptive CFL.

## Further Reading

- [`cfd-core` source](../../crates/cfd-core/src/)
- [`cfd-2d` source](../../crates/cfd-2d/src/)
- [`cfd-validation` source](../../crates/cfd-validation/src/)
- [Atlas Dependency Map](appendix_dependencies.md)
- Ghia, Ghia & Shin, "High-Re solutions for incompressible flow using the
  Navier-Stokes equations and a multigrid method", *J. Comput. Phys.* 48,
  387–411 (1982).
- White, *Viscous Fluid Flow*, 3rd ed., McGraw-Hill (2006) — Chapter 3
  for nondimensionalization.
