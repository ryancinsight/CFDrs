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



/// 2D CFD domain-specific prelude for advanced grid operations
///
/// This prelude exports 2D-specific functionality not available in the main prelude.
/// Use this when working extensively with 2D grids and need access to specialized
/// grid operations, advanced solver configurations, or detailed boundary handling.
///
/// For basic 2D functionality, prefer `cfd_suite::prelude::*`.
pub mod prelude {
    // === Advanced Grid Operations ===
    // Detailed grid functionality beyond basic StructuredGrid2D
    pub use crate::grid::{GridEdge, GridIterator};

    // === Specialized Solver Components ===
    // Advanced solver configurations and specialized methods
    pub use crate::{
        fdm::{AdvectionDiffusionSolver},
        fvm::{FluxScheme, Face},
        lbm::{D2Q9},
    };
}