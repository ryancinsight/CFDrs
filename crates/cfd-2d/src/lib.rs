//! 2D CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod grid;
pub mod schemes;
pub mod pressure_velocity;

// TODO: Implement these exports
// pub use fdm::{FiniteDifferenceMethod, FdmSolver};
// pub use fvm::{FiniteVolumeMethod, FvmSolver};
// pub use lbm::{LatticeBoltzmannMethod, LbmSolver, D2Q9};
// pub use grid::{Grid2D, StructuredGrid2D, UnstructuredGrid2D};
// pub use schemes::{AdvectionScheme, DiffusionScheme, Upwind, CentralDifference};
// pub use pressure_velocity::{SIMPLE, PISO, ProjectionMethod};

/// Common 2D CFD types and traits
pub mod prelude {
    // TODO: Add exports when implemented
    // pub use crate::{
    //     fdm::FdmSolver,
    //     fvm::FvmSolver,
    //     lbm::LbmSolver,
    //     grid::Grid2D,
    //     schemes::{Upwind, CentralDifference},
    //     pressure_velocity::SIMPLE,
    // };
}