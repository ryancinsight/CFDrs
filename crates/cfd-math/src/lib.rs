//! # CFD-Math: Advanced Mathematical Methods for Computational Fluid Dynamics
//!
//! This crate provides high-performance, numerically robust implementations of advanced
//! mathematical methods essential for computational fluid dynamics simulations.
//!
//! ## Core Mathematical Frameworks
//!
//! ### Linear Algebra & Sparse Matrices
//! - **CSR/CSC Storage**: Memory-efficient sparse matrix representations
//! - **Iterative Solvers**: GMRES, BiCGSTAB, Conjugate Gradient with preconditioning
//! - **Preconditioners**: ILU(k), AMG, Jacobi, SOR for solver acceleration
//!
//! ### High-Order Methods
//! - **Spectral Methods**: Fourier/Chebyshev expansions for DNS accuracy
//! - **Discontinuous Galerkin**: High-order finite elements with limiting
//! - **WENO Schemes**: Shock-capturing with 5th-order accuracy in smooth regions
//!
//! ### Numerical Analysis
//! - **Differentiation**: Finite difference, spectral differentiation
//! - **Integration**: Quadrature rules, composite integration
//! - **Interpolation**: Lagrange, cubic splines, ENO reconstruction
//!
//! ### Vectorization & Performance
//! - **SIMD Operations**: AVX2/SSE4.1 optimized kernels
//! - **Parallel Algorithms**: Rayon-based parallelism
//! - **Cache-Optimized**: Block-based algorithms for modern architectures
//!
//! ## Fundamental Mathematical Theorems
//!
//! ### Linear System Theory
//!
//! **Theorem (Convergence of Krylov Methods)**: For symmetric positive definite matrices,
//! the conjugate gradient method converges in at most n steps in exact arithmetic.
//!
//! **Theorem (GMRES Convergence)**: The GMRES method finds the best approximation
//! in the Krylov subspace K_m(A, r₀), providing optimal convergence properties.
//!
//! ### Preconditioning Theory
//!
//! **Theorem (ILU Stability)**: Incomplete LU factorization provides a stable
//! preconditioner when the fill-in level k is chosen appropriately for the matrix bandwidth.
//!
//! **Theorem (AMG Convergence)**: Algebraic multigrid methods achieve O(1) convergence
//! independent of problem size for elliptic operators.
//!
//! ### High-Order Accuracy
//!
//! **Theorem (Spectral Convergence)**: For smooth solutions, spectral methods achieve
//! exponential convergence O(e^(-c N)) where N is the number of modes.
//!
//! **Theorem (WENO Accuracy)**: WENO5 schemes achieve 5th-order accuracy in smooth regions
//! and maintain 3rd-order accuracy near discontinuities.
//!
//! **Theorem (DG Stability)**: Discontinuous Galerkin methods with appropriate limiting
//! are stable and high-order accurate for nonlinear conservation laws.
//!
//! ### Vectorization Theory
//!
//! **Theorem (SIMD Efficiency)**: SIMD operations provide theoretical speedup of up to
//! the vector width (4x for SSE4.1, 8x for AVX2) for appropriate data layouts.
//!
//! ## Implementation Philosophy
//!
//! ### Numerical Stability
//! - **Mixed Precision**: f64 for accuracy, f32 for performance where appropriate
//! - **Condition Number Control**: Preconditioning to maintain solver stability
//! - **Error Accumulation Control**: Careful ordering to minimize floating-point errors
//!
//! ### Performance Optimization
//! - **Memory Layout**: Cache-friendly data structures and access patterns
//! - **Algorithmic Complexity**: Optimal O(n) or O(n log n) scaling where possible
//! - **Parallel Scalability**: Thread-safe implementations with minimal synchronization
//!
//! ### Scientific Correctness
//! - **Conservation Properties**: Exact conservation of mass, momentum, energy
//! - **Symmetry Preservation**: Maintenance of physical symmetries in discretization
//! - **Error Estimation**: Built-in error indicators and adaptive refinement support
//!
//! ## Advanced Features
//!
//! ### Adaptive Methods
//! - **h-Adaptivity**: Mesh refinement based on error estimates
//! - **p-Adaptivity**: Polynomial order adaptation in spectral/DG methods
//! - **hp-Adaptivity**: Combined refinement strategies
//!
//! ### Specialized Solvers
//! - **Multigrid Methods**: Geometric and algebraic multigrid hierarchies
//! - **Matrix-Free Methods**: Reduced memory footprint for large problems
//! - **Domain Decomposition**: Parallel solution of large-scale problems
//!
//! ### Validation & Testing
//! - **Method of Manufactured Solutions**: Exact error verification
//! - **Convergence Studies**: Systematic accuracy assessment
//! - **Benchmark Problems**: Standardized test cases with known solutions
//!
//! ## References
//!
//! - **Iterative Methods**: Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*
//! - **Preconditioning**: Benzi, M. (2002). *Preconditioning Techniques for Large Linear Systems*
//! - **Spectral Methods**: Boyd, J.P. (2001). *Chebyshev and Fourier Spectral Methods*
//! - **DG Methods**: Cockburn, B. & Shu, C.-W. (2001). *Runge-Kutta Discontinuous Galerkin Methods*
//! - **WENO Methods**: Shu, C.-W. (1999). *High Order WENO Schemes for Convection Dominated Problems*
//! - **SIMD Programming**: Fog, A. (2012). *Optimizing software in C++*

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::needless_range_loop)] // Explicit indexing clearer for numerical algorithms
#![allow(clippy::too_many_lines)] // Complex numerical algorithms need detailed implementation
// CFD numerical computation allows
#![allow(clippy::similar_names)] // Mathematical variables often have similar names (x,y,z; i,j,k)
#![allow(clippy::cast_precision_loss)] // Precision loss acceptable for performance in numerical code
#![allow(clippy::cast_possible_truncation)] // Array indices and loop counters are typically small
#![allow(clippy::unused_self)] // Trait methods maintain interface consistency
#![allow(clippy::must_use_candidate)] // Mathematical utilities often used in larger expressions
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
#![allow(clippy::used_underscore_binding)] // Underscore prefixed bindings used for intentional partial use

pub mod differentiation;
pub mod error;
pub mod high_order;
pub mod integration;
pub mod interpolation;
pub mod iterators;
pub mod linear_solver;
pub mod performance_monitor;
pub mod pressure_velocity;
pub mod simd;
pub mod sparse;
pub mod time_stepping;

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
