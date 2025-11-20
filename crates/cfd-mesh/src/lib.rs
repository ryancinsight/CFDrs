//! Mesh handling and generation.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// Mesh generation and geometry allows
#![allow(clippy::similar_names)] // Geometric variables (x,y,z; i,j,k) naturally similar
#![allow(clippy::cast_precision_loss)] // Acceptable in geometric calculations
#![allow(clippy::must_use_candidate)] // Geometric utilities often used in expressions
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::ptr_arg)] // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations

pub mod connectivity;
#[cfg(feature = "csg")]
pub mod csg;
pub mod error;
pub mod geometry;
pub mod grid;
pub mod mesh;
pub mod quality;
pub mod refinement;
pub mod topology;

// The public modules are the primary API.
// Users should access types through the module hierarchy:
//   use cfd_mesh::mesh::Mesh;
//   use cfd_mesh::grid::StructuredGridBuilder;
//   use cfd_mesh::quality::aspect_ratio;
// The prelude provides convenient shortcuts for common types.

/// Common mesh types and traits
pub mod prelude {
    #[cfg(feature = "csg")]
    pub use crate::csg::CsgError;
    pub use crate::{
        connectivity::Connectivity,
        geometry::Geometry,
        mesh::{Mesh, MeshStatistics},
        topology::{Cell, Edge, Face, Vertex},
    };
}
