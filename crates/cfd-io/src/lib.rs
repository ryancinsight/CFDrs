//! Input/output functionality for CFD simulations.
//!
//! This crate provides file format support for reading and writing
//! CFD simulation data.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// I/O and file format allows
#![allow(clippy::similar_names)] // File format variables often similar
#![allow(clippy::must_use_candidate)] // I/O utilities typically return values that are used
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

pub mod binary;
pub mod checkpoint;
pub mod csv;
pub mod vtk;

// The API is now the public module hierarchy. This provides a clear,
// self-documenting structure for users.
// Example usage:
//   use cfd_io::vtk::VtkWriter;
//   use cfd_io::checkpoint::CheckpointManager;
//
// This hierarchical structure prevents namespace pollution and makes
// it clear which file format each type belongs to.
