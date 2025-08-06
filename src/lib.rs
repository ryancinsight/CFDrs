//! Computational Fluid Dynamics Simulation Suite
//!
//! A comprehensive CFD simulation framework in Rust supporting 1D, 2D, and 3D simulations
//! with a plugin-based architecture.

#![warn(missing_docs)]

// Re-export all sub-crates
pub use cfd_core as core;
pub use cfd_math as math;
pub use cfd_io as io;
pub use cfd_mesh as mesh;
pub use cfd_1d as d1;
pub use cfd_2d as d2;
pub use cfd_3d as d3;
pub use cfd_validation as validation;

/// Prelude module for convenient imports
pub mod prelude {
    pub use cfd_core::prelude::*;
    pub use cfd_math::prelude::*;
    pub use cfd_io::prelude::*;
    pub use cfd_mesh::prelude::*;
    pub use cfd_1d::prelude::*;
    pub use cfd_2d::prelude::*;
    pub use cfd_3d::prelude::*;
    pub use cfd_validation::prelude::*;
}