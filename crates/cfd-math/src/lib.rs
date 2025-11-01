//! Mathematical utilities and numerical methods for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::needless_range_loop)]      // Explicit indexing clearer for numerical algorithms
#![allow(clippy::too_many_lines)]           // Complex numerical algorithms need detailed implementation
// CFD numerical computation allows
#![allow(clippy::similar_names)]           // Mathematical variables often have similar names (x,y,z; i,j,k)
#![allow(clippy::cast_precision_loss)]     // Precision loss acceptable for performance in numerical code
#![allow(clippy::cast_possible_truncation)] // Array indices and loop counters are typically small
#![allow(clippy::unused_self)]             // Trait methods maintain interface consistency
#![allow(clippy::must_use_candidate)]      // Mathematical utilities often used in larger expressions
#![allow(clippy::missing_errors_doc)]      // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)]      // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)]          // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)]      // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)]      // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)]               // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)]        // Result types maintained for API consistency
#![allow(clippy::items_after_statements)]  // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)]       // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)]      // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)]            // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)]  // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)]  // Builder patterns used internally
#![allow(clippy::ptr_arg)]                 // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)]  // CFD-specific trait implementations
// CFD numerical computation allows
#![allow(clippy::similar_names)]           // Mathematical variables often have similar names (x,y,z; i,j,k)
#![allow(clippy::cast_precision_loss)]     // Precision loss acceptable for performance in numerical code
#![allow(clippy::cast_possible_truncation)] // Array indices and loop counters are typically small
#![allow(clippy::unused_self)]             // Trait methods maintain interface consistency
#![allow(clippy::must_use_candidate)]      // Mathematical utilities often used in larger expressions
#![allow(clippy::missing_errors_doc)]      // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)]      // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)]          // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)]      // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)]      // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)]               // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)]        // Result types maintained for API consistency
#![allow(clippy::items_after_statements)]  // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)]       // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)]      // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)]            // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)]  // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)]  // Builder patterns used internally
#![allow(clippy::ptr_arg)]                 // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)]  // CFD-specific trait implementations
#![allow(clippy::used_underscore_binding)]  // Underscore prefixed bindings used for intentional partial use

pub mod differentiation;
pub mod error;
// pub mod high_order; // Temporarily disabled - needs implementation
pub mod integration;
pub mod interpolation;
pub mod iterators;
pub mod linear_solver;
pub mod preconditioners;
pub mod simd;
pub mod sparse;
pub mod vector_ops;
pub mod vectorization;
pub mod cfd_simd;

// --- Curated Top-Level API ---
// Only expose a very small number of absolutely fundamental traits or structs.
// Users should interact with the module hierarchy for most types.

pub use self::interpolation::Interpolation;
pub use self::sparse::SparseMatrix;

// The primary API is through the public modules themselves.
// This creates a hierarchical, self-documenting structure.
// Example usage:
//   use cfd_math::linear_solver::BiCGSTAB;
//   use cfd_math::interpolation::CubicSplineInterpolation;
//   use cfd_math::integration::GaussQuadrature;

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        differentiation::FiniteDifference,
        integration::Quadrature,
        interpolation::{Interpolation, LinearInterpolation},
        linear_solver::ConjugateGradient,
        sparse::{SparseMatrix, SparseMatrixBuilder},
    };
}
