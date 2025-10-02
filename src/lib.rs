//! Computational Fluid Dynamics Simulation Suite
//!
//! A comprehensive CFD simulation framework in Rust supporting 1D, 2D, and 3D simulations
//! with a plugin-based architecture.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// CFD Suite unified configuration - aligning with all workspace crates
#![allow(clippy::similar_names)]           // CFD variables across dimensions
#![allow(clippy::cast_precision_loss)]     // Numerical computing requires performance trade-offs
#![allow(clippy::cast_possible_truncation)] // Grid indices typically small
#![allow(clippy::unused_self)]             // Trait method interface consistency
#![allow(clippy::must_use_candidate)]      // Computational utilities context-dependent
#![allow(clippy::missing_errors_doc)]      // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)]      // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)]          // CFD indexing patterns
#![allow(clippy::cast_possible_wrap)]      // Grid index wrap-around acceptable
#![allow(clippy::too_many_arguments)]      // CFD functions need many physical parameters
#![allow(clippy::float_cmp)]               // Necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)]        // API consistency for Result types
#![allow(clippy::items_after_statements)]  // Readability in numerical code
#![allow(clippy::many_single_char_names)]       // Mathematical notation standard (i,j,k,x,y,z)
#![allow(clippy::unreadable_literal)]      // Physical constants require precision
#![allow(clippy::redundant_closure_for_method_calls)] // Pipeline readability
#![allow(clippy::doc_markdown)]            // Math notation formatting
#![allow(clippy::needless_pass_by_value)]  // Copy types idiom
#![allow(clippy::return_self_not_must_use)]  // Builder pattern support
#![allow(clippy::ptr_arg)]                 // API compatibility
#![allow(clippy::should_implement_trait)]  // CFD-specific implementations

// Re-export all sub-crates
pub use cfd_1d as d1;
pub use cfd_2d as d2;
pub use cfd_3d as d3;
pub use cfd_core as core;
pub use cfd_io as io;
pub use cfd_math as math;
pub use cfd_mesh as mesh;
pub use cfd_validation as validation;

pub mod compute_unified;

/// Unified prelude module - Single Source of Truth for all common CFD functionality
///
/// This module provides a comprehensive set of imports for CFD simulations across all dimensions.
/// It follows the SSOT (Single Source of Truth) principle by being the primary interface for users.
///
/// # Usage
/// ```rust
/// use cfd_suite::prelude::*;
/// use cfd_suite::core::error::Result;
/// use cfd_suite::d1::Network;
/// use cfd_suite::d2::grid::StructuredGrid2D;
///
/// fn main() -> Result<()> {
///     // Create fluid properties using proper API
///     let fluid = ConstantPropertyFluid::<f64>::water_20c()?;
///     // Create network graph first, then network with fluid
///     let graph = cfd_suite::d1::network::NetworkGraph::new();
///     let mut network = Network::new(graph, fluid.clone());
///     let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0);
///     Ok(())
/// }
/// ```
pub mod prelude {
    //! Common imports for CFD simulations

    // Core functionality
    pub use cfd_core::prelude::*;

    // Math operations
    pub use cfd_math::prelude::*;

    // Mesh operations - exclude Edge to avoid conflict
    pub use cfd_mesh::prelude::{Cell, Connectivity, Face, Geometry, Mesh, Vertex};
    // Re-export mesh Edge with qualifier to avoid ambiguity
    pub use cfd_mesh::topology::Edge as MeshEdge;

    // I/O operations - use proper module paths
    pub use cfd_io::vtk::{VtkMesh, VtkMeshBuilder, VtkReader, VtkWriter};

    // 1D solver exports - exclude Edge to avoid conflict
    pub use cfd_1d::prelude::{
        Channel, ChannelGeometry, ChannelType, ComponentFactory, CrossSection, DarcyWeisbachModel,
        EdgeProperties, EdgeType, FlowConditions, FlowRegime, FlowSensor, FlowState, Micromixer,
        MixerType, NetworkAnalysisResult, Node, NodeProperties, NodeType, PerformanceMetrics,
        PumpType, ResistanceCalculator, ResistanceModelFactory, SensorType, SurfaceProperties,
        ValveType, Wettability,
    };
    // Re-export 1D Edge with qualifier
    pub use cfd_1d::network::Edge as NetworkEdge;

    // 3D CFD functionality
    pub use cfd_3d::*;

    // Note: cfd_2d and cfd_3d don't have prelude modules yet

    // Validation tools - use proper module paths
    pub use cfd_validation::convergence::{
        ConvergenceStudy, GridConvergenceIndex, RichardsonExtrapolation,
    };
    pub use cfd_validation::error_metrics::{ErrorAnalysis, ErrorMetric};
}
