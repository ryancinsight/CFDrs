//! # cfd-2d — Full 2D Navier-Stokes PDE Solver
//!
//! `cfd-2d` solves the **incompressible 2D Navier-Stokes equations** on
//! structured and unstructured grids. Unlike `cfd-1d`, which reduces each
//! channel to a single lumped resistance, `cfd-2d` resolves the full spatial
//! velocity and pressure fields (u, v, p) at every grid cell.
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
//! | Solver        | Module                    | Method                                   |
//! |---------------|---------------------------|------------------------------------------|
//! | FDM           | `solvers::fdm`            | Finite Difference (Poisson, diffusion)   |
//! | FVM           | `solvers::fvm`            | Finite Volume (flux-based)               |
//! | LBM           | `solvers::lbm`            | Lattice-Boltzmann D2Q9                   |
//! | SIMPLE        | `solvers::simple`         | Semi-Implicit pressure-velocity coupling |
//! | PISO          | `piso_algorithm`          | Pressure-Implicit Split Operator         |
//! | SIMPLEC/PIMPLE| `simplec_pimple`          | Extended SIMPLE variants                 |
//! | Poiseuille    | `solvers::poiseuille`     | Analytical Poiseuille + non-Newtonian    |
//! | Bifurcation   | `solvers::bifurcation_flow` | 2D bifurcating channel flow            |
//! | Serpentine    | `solvers::serpentine_flow`  | 2D serpentine channel flow             |
//! | Venturi       | `solvers::venturi_flow`     | 2D Venturi constriction flow           |
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
//!
//! ## Relationship with `cfd-schematics`
//!
//! `cfd-2d` does **not** currently depend on `cfd-schematics`. The 2D
//! simulation domain is a continuous PDE field on a grid — it is not a lumped
//! network graph. The "2D" in `cfd-schematics` refers to the **layout plane**
//! of a chip schematic (x, y in mm), not a PDE simulation domain.
//!
//! `cfd-schematics` *can* support `cfd-2d` workflows in two ways:
//! 1. **Geometry seeding**: `ChannelSystem` centrelines and bounding boxes can
//!    define the domain extents and inlet/outlet positions for a 2D grid.
//! 2. **Result overlay**: `AnalysisOverlay` can project cross-section-averaged
//!    2D solver results back onto the schematic for design-level inspection.
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


#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// 2D CFD simulation allows
#![allow(clippy::similar_names)] // CFD variables (u,v,p; nx,ny; dx,dy; i,j) often similar
#![allow(clippy::cast_precision_loss)] // Performance-critical numerical loops
#![allow(clippy::cast_possible_truncation)] // Grid indices and array sizes typically small
#![allow(clippy::unused_self)] // Solver trait methods maintain consistent interfaces
#![allow(clippy::must_use_candidate)] // Solver utilities and getters used in computational contexts
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::ptr_arg)] // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations
#![allow(clippy::too_many_lines)] // Complex solver implementations require detailed methods
#![allow(clippy::needless_range_loop)] // Explicit indexing clearer for multi-dimensional CFD arrays
#![allow(clippy::struct_field_names)] // Field names like field_* common in computational contexts
#![allow(clippy::used_underscore_binding)] // Underscore prefixed bindings used for intentional partial use

// Core modules
pub mod constants;
pub mod error;
pub mod fields;
pub mod grid;
pub mod problem;

// Domain-organized modules
pub mod discretization;
pub mod physics;
pub mod solvers;

// Algorithm modules
pub mod piso_algorithm;
pub mod pressure_velocity;
pub mod simplec_pimple;

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
