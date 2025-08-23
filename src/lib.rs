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

/// Unified prelude module - Single Source of Truth for all common CFD functionality
///
/// This module provides a comprehensive set of imports for CFD simulations across all dimensions.
/// It follows the SSOT (Single Source of Truth) principle by being the primary interface for users.
///
/// # Usage
/// ```rust
/// use cfd_suite::prelude::*;
/// use cfd_suite::core::Result;
/// use cfd_suite::d1::Network;
/// use cfd_suite::d2::StructuredGrid2D;
///
/// fn main() -> Result<()> {
///     // Now you have access to all commonly used CFD functionality
///     let fluid = Fluid::<f64>::water()?;
///     let mut network = Network::new(fluid);
///     let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
///     Ok(())
/// }
/// ```
pub mod prelude {
    //! Common imports for CFD simulations
    
    // Core types
    pub use cfd_core::prelude::*;
    
    // Math operations
    pub use cfd_math::prelude::*;
    
    // Mesh operations
    pub use cfd_mesh::prelude::*;
    
    // I/O operations - only export what's available
    pub use cfd_io::{VtkWriter, VtkReader, VtkMesh, VtkMeshBuilder};
    
    // Solver exports
    pub use cfd_1d::prelude::*;
    // Note: cfd_2d and cfd_3d don't have prelude modules yet
    // TODO: Add prelude modules to these crates for consistency
    
    // Validation tools
    pub use cfd_validation::prelude::*;
}