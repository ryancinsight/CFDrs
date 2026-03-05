//! Mesh smoothing algorithms for 2-D constrained Delaunay triangulations.
//!
//! ## Algorithms
//!
//! | Smoother | Guarantee | Use Case |
//! |----------|-----------|----------|
//! | [`LaplacianSmoother`] | Boundary-preserving convergence | General smoothing |
//! | [`AngleBasedSmoother`] | Monotone min-angle improvement | Quality-critical meshes |
//!
//! ## Usage
//!
//! Both smoothers operate on a `&mut Cdt` in place.  They preserve:
//! - The **triangulation connectivity** (no flips or vertex insertions).
//! - **Boundary vertex positions** (when `preserve_boundary = true`).
//!
//! After smoothing, the triangulation may no longer satisfy the strict
//! Delaunay empty-circumcircle property.  Re-running Lawson flips (via
//! `Cdt::from_pslg`) is recommended if strict Delaunay compliance is needed.
//!
//! ## Integration with `RuppertRefiner`
//!
//! See [`crate::application::delaunay::refinement::ruppert::RuppertRefiner`]
//! for the optional `post_smooth` configuration field that calls
//! `LaplacianSmoother::smooth` after refinement completes.

pub mod angle_based;
pub mod laplacian;

pub use angle_based::AngleBasedSmoother;
pub use laplacian::LaplacianSmoother;
