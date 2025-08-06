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
pub mod simple;

// Export implemented functionality
pub use grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge, GridIterator};
pub use fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig};
pub use fvm::{FvmSolver, FvmConfig, FluxScheme};
pub use lbm::{LbmSolver, LbmConfig, D2Q9};
pub use simple::{SimpleSolver, SimpleConfig};

// TODO: Implement these exports
// pub use schemes::{AdvectionScheme, DiffusionScheme, Upwind, CentralDifference};

/// Common 2D CFD types and traits
pub mod prelude {
    // Export implemented functionality
    pub use crate::{
        grid::{Grid2D, StructuredGrid2D, BoundaryType, GridEdge},
        fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig},
        fvm::{FvmSolver, FvmConfig, FluxScheme},
        lbm::{LbmSolver, LbmConfig, D2Q9},
        simple::{SimpleSolver, SimpleConfig},
    };

    // TODO: Add exports when implemented
    // pub use crate::{
    //     schemes::{Upwind, CentralDifference},
    // };
}