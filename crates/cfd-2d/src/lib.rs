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

// Export implemented functionality
pub use grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge, GridIterator};
pub use fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig};

// TODO: Implement these exports
// pub use fvm::{FiniteVolumeMethod, FvmSolver};
// pub use lbm::{LatticeBoltzmannMethod, LbmSolver, D2Q9};
// pub use schemes::{AdvectionScheme, DiffusionScheme, Upwind, CentralDifference};
// pub use pressure_velocity::{SIMPLE, PISO, ProjectionMethod};

/// Common 2D CFD types and traits
pub mod prelude {
    // Export implemented functionality
    pub use crate::{
        grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge},
        fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig},
    };

    // TODO: Add exports when implemented
    // pub use crate::{
    //     fvm::FvmSolver,
    //     lbm::LbmSolver,
    //     schemes::{Upwind, CentralDifference},
    //     pressure_velocity::SIMPLE,
    // };
}