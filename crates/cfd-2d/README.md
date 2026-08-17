# cfd-2d

2D incompressible Navier-Stokes solvers for the
[CFDrs](https://github.com/ryancinsight/CFDrs) simulation suite.

Where `cfd-1d` reduces each channel to a lumped resistance, `cfd-2d` resolves the
velocity and pressure fields `(u, v, p)` at every grid cell.

## Solvers

| Family | Module | Method |
|---|---|---|
| FDM | `solvers::fdm` | Finite difference (Poisson, diffusion) |
| FVM | `solvers::fvm` | Finite volume (primary production family) |
| LBM | `solvers::lbm` | Lattice-Boltzmann D2Q9 cross-check path |
| SIMPLE | `solvers::simple` | Semi-implicit pressure-velocity coupling |
| PISO | `piso_algorithm` | Pressure-implicit split operator |
| SIMPLEC / PIMPLE | `simplec_pimple` | Extended SIMPLE variants |
| Analytical | `solvers::poiseuille` | Poiseuille, including non-Newtonian |
| Device flows | `solvers::{bifurcation_flow, serpentine_flow, venturi_flow}` | Channel geometries |

Turbulence closures live in `physics::turbulence`: standard k-ε, k-ω SST (Menter),
detached eddy simulation, and a Reynolds stress model.

## Also here

- **`grid`** — `StructuredGrid2D`, `UnstructuredGrid2D`, boundary types, adaptive refinement.
- **`fields`** — `SimulationFields` containers for u, v, p, T.
- **`discretization`, `schemes`** — convection schemes and numerical fluxes.
- **`network`** — projects a `cfd-schematics` blueprint into per-channel 2D solves,
  using a 1D reference solve to configure boundary conditions.
- **`stability`** — CFL and stability analysis.

## Example

The crate-level Rustdoc carries a runnable quick-start doctest:

```bash
cargo test --doc -p cfd-2d
```

## Features

| Feature | Default | Effect |
|---|---|---|
| `gpu` | no | Enables `cfd-core/gpu` for accelerated kernels |

## Documentation

- [CFDrs book](https://ryancinsight.github.io/CFDrs/)
- [Workspace README](https://github.com/ryancinsight/CFDrs#readme)

## License

MIT OR Apache-2.0
