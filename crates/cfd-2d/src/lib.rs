//! 2D CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod energy;
pub mod fdm;
pub mod field;  // Current field abstraction module
pub mod fields;
pub mod fvm;
pub mod grid;
pub mod lbm;
pub mod momentum;
pub mod convection;
pub mod pressure_velocity;  // Pressure-velocity coupling algorithm (STANDARD)
pub mod piso;

pub mod problem;
pub mod turbulence;
pub mod vorticity_stream;
pub mod schemes;

// Export implemented functionality
pub use grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge};
pub use energy::EnergyEquationSolver;
pub use turbulence::{KEpsilonModel, WallFunction};
pub use fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig};
pub use fvm::{FvmSolver, FvmConfig, FluxScheme};
pub use lbm::{LbmSolver, LbmConfig, D2Q9};
pub use pressure_velocity_coupling::{PressureVelocityCouplerSolver, PressureVelocityCouplingConfig};
pub use fields::{Field2D, SimulationFields};
pub use momentum::{MomentumSolver, MomentumComponent, MomentumCoefficients};
pub use problem::{IncompressibleFlowProblem, IncompressibleFlowSolution};
pub use convection::{ConvectionScheme, ConvectionSchemeFactory, FirstOrderUpwind, CentralDifference, HybridScheme, PowerLawScheme, QuickScheme};
pub use schemes::{SpatialScheme, FluxLimiter, TimeScheme, FiniteDifference};



// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface