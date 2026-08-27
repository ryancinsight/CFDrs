//! # cfd-2d — Full 2D Navier-Stokes PDE Solver
//!
//! `cfd-2d` solves the **incompressible 2D Navier-Stokes equations** for
//! laminar, low-Mach millifluidic flows on structured grids. Unlike `cfd-1d`,
//! which reduces each channel to a single lumped resistance, `cfd-2d`
//! resolves the spatial velocity and pressure fields `(u, v, p)` at every
//! grid cell and uses schematics-driven projection to map blueprint channels
//! into solver masks and boundary-aware metrics.
//!
//! ## Physical Model
//!
//! The governing equations are the incompressible Navier-Stokes equations:
//!
//! ```text
//!   ∂u/∂t + (u·∇)u = -∇p/ρ + ν·∇²u + f    (momentum)
//!   ∇·u = 0                                  (continuity)
//! ```
//!
//! These are solved on a 2D spatial domain (x, y) with appropriate boundary
//! conditions (no-slip walls, inlet velocity, outlet pressure, etc.).
//!
//! ## Solvers
//!
//! | Solver        | Module                      | Method                                     |
//! |---------------|-----------------------------|--------------------------------------------|
//! | FDM           | `solvers::fdm`              | Finite Difference (Poisson, diffusion)     |
//! | FVM           | `solvers::fvm`              | Finite Volume (primary production family)  |
//! | LBM           | `solvers::lbm`              | Lattice-Boltzmann D2Q9 cross-check path     |
//! | SIMPLE        | `solvers::simple`           | Semi-Implicit pressure-velocity coupling    |
//! | PISO          | `piso_algorithm`            | Pressure-Implicit Split Operator            |
//! | SIMPLEC/PIMPLE| `simplec_pimple`            | Extended SIMPLE variants                    |
//! | Poiseuille    | `solvers::poiseuille`       | Analytical Poiseuille + non-Newtonian       |
//! | Bifurcation   | `solvers::bifurcation_flow`  | 2D bifurcating channel flow                 |
//! | Serpentine    | `solvers::serpentine_flow`   | 2D serpentine channel flow                  |
//! | Venturi       | `solvers::venturi_flow`      | 2D Venturi constriction flow                |
//!
//! ## Turbulence Models
//!
//! - `physics::turbulence::k_epsilon` — Standard k-ε
//! - `physics::turbulence::k_omega_sst` — k-ω SST (Menter)
//! - `physics::turbulence::des` — Detached Eddy Simulation
//! - `physics::turbulence::reynolds_stress` — Reynolds Stress Model
//!
//! ## Relationship with `cfd-1d`
//!
//! | Aspect              | `cfd-1d`                          | `cfd-2d`                              |
//! |---------------------|-----------------------------------|---------------------------------------|
//! | Governing equations | Kirchhoff / Hagen-Poiseuille      | Navier-Stokes (u, v, p fields)        |
//! | Spatial resolution  | Per-channel scalar (Q, ΔP)        | Per-cell vector/scalar field          |
//! | Grid                | Graph (nodes + edges)             | Structured / unstructured 2D grid     |
//! | Turbulence          | Not applicable (Re ≪ 1)           | k-ε, k-ω SST, DES, Reynolds stress   |
//! | Cost                | Microseconds                      | Seconds to hours                      |
//! | Use case            | Network design, flow distribution | Detailed velocity/pressure fields     |
//!
//! A typical workflow is: run `cfd-1d` first to find the pressure distribution
//! across the network, then use those pressures as inlet/outlet boundary
//! conditions for a `cfd-2d` simulation of a single critical channel segment.
//! The `network` module now follows that pattern directly for
//! `cfd-schematics::NetworkBlueprint` inputs by building a normalized 1D
//! reference solve, projecting each channel path into a 2D mask, and attaching
//! per-node, per-channel, and projection summaries to the channel solves.
//!
//! ## Relationship with `cfd-schematics`
//!
//! `cfd-2d` consumes `cfd-schematics::NetworkBlueprint` in its `network`
//! module. The blueprint remains the topology and geometry single source of
//! truth, while `cfd-2d` builds one PDE solve per channel, rasterizes routed
//! paths into solver masks, and retains the 1D trace used to configure those
//! solves. `solve_projected` returns the compatibility result together with the
//! projection metadata for audits and benchmarks.
//!
//! # Modules
//! - **grid**: `StructuredGrid2D`, `UnstructuredGrid2D`, boundary types, refinement
//! - **fields**: `SimulationFields` (u, v, p, T scalar/vector field containers)
//! - **problem**: `IncompressibleFlowProblem`, `IncompressibleFlowSolution`
//! - **solvers**: FDM, FVM, LBM, SIMPLE, Poiseuille, bifurcation, serpentine, venturi
//! - **physics**: Momentum, energy, turbulence, immersed boundary, vorticity-stream
//! - **discretization**: Convection schemes, extended stencils
//! - **piso_algorithm**: PISO predictor/corrector loop
//! - **pressure_velocity**: Pressure-velocity coupling coefficients
//! - **simplec_pimple**: SIMPLEC and PIMPLE algorithm variants
//! - **schemes**: Numerical flux schemes
//! - **stability**: CFL and stability analysis utilities
//! - **constants**: Physical and numerical constants
//!
//! # Quick Start
//!
//! Construct a validated fluid and a structured grid, and derive a physical
//! quantity from them. Both constructors are fallible and share
//! [`cfd_core::error::Error`], so a single `?` chain covers them.
//!
//! ```
//! use cfd_2d::grid::{Grid2D, StructuredGrid2D};
//! use cfd_core::physics::fluid::{ConstantFluid, ConstantPropertyFluid};
//!
//! # fn main() -> Result<(), cfd_core::error::Error> {
//! // Water at 20 °C. Properties are range-checked during construction.
//! let water = ConstantPropertyFluid::<f64>::water_20c()?;
//! let nu = water.dynamic_viscosity().into_base() / water.density().into_base();
//!
//! // 64x64 structured grid over the unit square [0,1] x [0,1].
//! let grid = StructuredGrid2D::<f64>::unit_square(64, 64)?;
//! assert_eq!(grid.num_cells(), 64 * 64);
//!
//! // Kinematic viscosity of water at 20 °C: 1.002e-3 / 998.2 ~= 1.004e-6 m^2/s.
//! assert!((nu - 1.004e-6).abs() < 1e-9);
//! # Ok(())
//! # }
//! ```
//!
//! This covers construction and property access only. Driving a solver to
//! convergence requires a boundary-condition set and a time or iteration loop;
//! see the [`solvers`] module and the repository `examples/` directory.
//!
//! # Discrete Invariants
//!
//! For the supported incompressible formulations, the implementation is expected to preserve
//! three crate-level invariants on every converged update:
//!
//! 1. Discrete mass conservation after the pressure-correction step.
//! 2. Non-negative modeled transport coefficients such as eddy viscosity and scalar diffusivity.
//! 3. Bounded face reconstruction when TVD or other limited high-resolution schemes are selected
//!    within their stated stability envelope.
//!
//! **Proof sketch**:
//! The SIMPLE, SIMPLEC, PISO, and PIMPLE paths assemble a momentum predictor and then solve a
//! pressure-correction problem whose flux correction enforces the discrete continuity constraint
//! on each control volume. Rhie-Chow interpolation removes collocated pressure checkerboarding so
//! the corrected flux field remains consistent with the pressure solve. Turbulence and transport
//! closures compute viscosity-like quantities from non-negative state variables, and the bounded
//! convection schemes limit reconstructed interface values so new extrema are not introduced when
//! their CFL assumptions hold.

#![warn(missing_docs)]
// 2D CFD simulation allows

// Core modules
pub mod constants;
pub mod fields;
pub mod grid;
pub mod problem;
pub(crate) mod scalar;

// Domain-organized modules
pub mod discretization;
pub mod physics;
pub mod solvers;

// Algorithm modules
pub mod piso_algorithm;
pub mod pressure_velocity;
pub mod simplec_pimple;

pub mod network;
pub mod schemes;
pub mod stability;

// The crate's public API is its module hierarchy.
// Users should access types with clear, logical paths:
//   use cfd_2d::solvers::fvm::FvmSolver;
//   use cfd_2d::physics::turbulence::KEpsilonModel;
//   use cfd_2d::discretization::ConvectionScheme;
//   use cfd_2d::fields::SimulationFields;
//   use cfd_2d::grid::StructuredGrid2D;
//
// This hierarchical structure is self-documenting and aligns with Rust best practices.

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface
