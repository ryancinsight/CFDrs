# CFDrs Architecture and Problem Setup

CFDrs is a workspace of CFD-focused crates organized so that *physics,
geometry, and solver* concerns live in **distinct** crates and compose through
workspace-level contracts.  This Part introduces the architecture and points
at the example scenarios that exercise every contract.

## Crate Map

| Crate | Role |
|---|---|
| `cfd-core` | Scalar seams, dimensional types, boundary conditions, error types |
| `cfd-math` | Linear algebra, sparse matrices, iterative solvers, FFT bridge |
| `cfd-1d` | 1-D blood/Reynolds models, microfluidic screens, hemolysis pipelines |
| `cfd-2d` | 2-D Stencil solvers, biomass CFD, CSG-coupled channels |
| `cfd-3d` | FEM and spectral 3-D solvers, multi-region flow |
| `cfd-validation` | Reference phantoms and orthogonal oracles (Ghia, Hagen–Poiseuille, ...) |
| `cfd-schematics` | CSG primitives and curve generators (venturi, bifurcation, serpentine) |
| `cfd-io` | HDF5 / NetCDF / CSV output adapters |
| `cfd-optim` | Optimization audit harness (Ray tracing milestone metas) |
| `cfd-python` | PyO3 bindings to Python-friendly data products |

Every solver crate depends on `cfd-core` and `cfd-math`; geometry crates
(`cfd-schematics`) are depended on by `cfd-{1,2,3}d`.  The cross-cutting
contracts — boundary conditions, time-stepping, state vectors — live in
`cfd-core` and are referenced by example chapters throughout this book.

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
Every part of this book references one of those three components.

## The CFDrs Simulation Lifecycle

```
+------------+   +-----------+   +----------+   +----------+
| Geometry   |-->| Mesh/Grid |-->| Boundary |-->| Initial  |
| (CSG)      |   |           |   | wiring   |   | state    |
+------------+   +-----------+   +----------+   +----------+
                                                    |
                                                    v
+------------+   +-----------+   +----------+
| Posterior  |<--| Solve     |<--| Time     |
| analysis   |   | (linear   |   | step     |
|            |   |  iter.)   |   |          |
+------------+   +-----------+   +----------+
```

- **Geometry → Mesh** is handled in `cfd-schematics` (CSG) and `cfd-{2,3}d/grid.rs`.
- **Boundary wiring** lives in [`cfd_core::boundary`].
- **Time stepping** is configurable per solver (`cfd-2d::time` /
  `cfd-3d::time` / adaptive in `cfd-2d::adaptive`).
- **Posterior** means DVH-equivalent (cavitation risk, shear envelope,
  Reynolds stress), reported in user-supplied hooks.

## Atlas Stack Interface

CFDrs is mid-migration to the Atlas stack.  Where a crate already declares
Atlas crates in `Cargo.toml`, the migration is considered **complete**; where
`ndarray` or `nalgebra` is still in the manifest, the crate is in
**bulk-migration** phase.  See [Atlas Dependency Map](appendix_dependencies.md)
for the per-crate table and [Atlas Stack Integration — Migration Reference]
for the migration chapters.  Migration progress also lives in
[`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md) and the per-crate roadmap
sections.

## Examples Referenced by This Chapter

- [`cfd_demo`](examples/cfd_demo.md) — minimum end-to-end setup.
- [`enhanced_cfd_demo`](examples/enhanced_cfd_demo.md) — diagnostics enabled.
- [`adaptive_time_stepping_demo`](examples/adaptive_time_stepping_demo.md)
  — demonstrates the time-stepping hook with adaptive CFL.

## Further Reading

- [`cfd-core` source](../../crates/cfd-core/src/)
- [`cfd-2d` source](../../crates/cfd-2d/src/)
- [Atlas Dependency Map](appendix_dependencies.md)
