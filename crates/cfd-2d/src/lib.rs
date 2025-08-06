//! 2D CFD solvers including FDM, FVM, and LBM implementations.

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

pub use fdm::{FiniteDifferenceMethod, FdmSolver};
pub use fvm::{FiniteVolumeMethod, FvmSolver};
pub use lbm::{LatticeBoltzmannMethod, LbmSolver, D2Q9};
pub use grid::{Grid2D, StructuredGrid2D, UnstructuredGrid2D};
pub use schemes::{AdvectionScheme, DiffusionScheme, Upwind, CentralDifference};
pub use pressure_velocity::{SIMPLE, PISO, ProjectionMethod};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        fdm::FdmSolver,
        fvm::FvmSolver,
        lbm::LbmSolver,
        grid::Grid2D,
        schemes::{Upwind, CentralDifference},
        pressure_velocity::SIMPLE,
    };
}