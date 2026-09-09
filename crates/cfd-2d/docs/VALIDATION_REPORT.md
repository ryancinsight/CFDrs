# cfd-2d Validation Report

## Solver Scope
- Primary production path: structured-grid finite-volume Navier-Stokes with SIMPLE, PISO, SIMPLEC, and PIMPLE variants.
- Secondary cross-check path: LBM D2Q9, kept for validation and comparison rather than as the default millifluidic workhorse.
- Applicable regime: laminar, low-Mach, low-Re millifluidic channels. Turbulence models remain available, but they are outside the default design envelope for the schematics-driven network path.

## Schematics Projection
- `cfd-schematics::NetworkBlueprint` remains the topology and geometry source of truth.
- The `network` module projects routed `ChannelSpec.path` data into 2D solver masks, boundary-aware extents, and per-channel projection summaries.
- Validation rejects unresolved routed crossings and degenerate channel paths before any 2D solve is started.
- The compatibility solve path remains per-channel; `solve_projected` returns the same solve result together with projection metadata for audits and benchmarks.

## Verified Invariants
- Anchored pressure-correction solves preserve discrete mass conservation to the linear-solver tolerance on the fluid domain.
- Rhie-Chow face interpolation suppresses checkerboard pressure modes on colocated grids when positive momentum diagonals and consistent stencils are used.
- TVD MUSCL reconstructions remain bounded only when the limiter stays in the admissible Sweby region and the explicit update satisfies the CFL assumptions.
- WENO reconstructions are not strictly TVD; boundedness relies on the nonlinear weights, the smoothness indicators, and a stable SSP time integrator.

## Test Coverage
- Blueprint validation tests cover routed-path admissibility, unresolved crossing rejection, and selective-tree blueprint acceptance.
- Network regression tests cover straight, serpentine, venturi, bifurcation, and selective-tree projection rasterization.
- Physics regression tests cover mass conservation, outlet-flow agreement, and wall-shear ordering on schematics-generated channels.
- Performance benchmarks now separate reference-trace construction, projection/build, core solve, and projected result wrapping.

## Research Basis
- Computational methods for inertial microfluidics: recent advances and future perspectives
  - https://www.nature.com/articles/s41378-025-00992-6
- Accelerated Computational Fluid Dynamics Simulations of Microfluidic Devices by Exploiting Higher Levels of Abstraction
  - https://www.cda.cit.tum.de/files/eda/2024_micromachines_mdpi_accelerated_computational_fluid_dynamics_simulations_of_microfluidic_devices_by_exploiting_higher_levels_of_abstraction.pdf
- A comprehensive uncertainty analysis of lattice Boltzmann flow simulations in a bifurcation geometry
  - https://www.nature.com/articles/s41598-025-30414-6.pdf
