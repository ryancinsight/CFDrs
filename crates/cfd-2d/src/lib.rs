//! 2D CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod constants;
pub mod fdm;
pub mod fvm;
pub mod grid;
pub mod lbm;
pub mod simple;
pub mod piso;
pub mod vorticity_stream;
pub mod schemes;

// Export implemented functionality
pub use grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge, GridIterator};
pub use fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig};
pub use fvm::{FvmSolver, FvmConfig, FluxScheme};
pub use lbm::{LbmSolver, LbmConfig, D2Q9};
pub use simple::{SimpleSolver, SimpleConfig};
pub use schemes::{SpatialScheme, FluxLimiter, TimeScheme, FiniteDifference};



// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface