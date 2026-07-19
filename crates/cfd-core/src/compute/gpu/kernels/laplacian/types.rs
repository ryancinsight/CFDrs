//! Types for the 2D Laplacian GPU kernel.

pub use hephaestus_wgpu::BoundaryCondition;

/// Deprecated alias for [`BoundaryCondition`].
///
/// Use [`BoundaryCondition`] directly. This alias is kept for backward
/// compatibility and will be removed in a future release.
#[deprecated(
    since = "0.18.0",
    note = "Use hephaestus_wgpu::BoundaryCondition (or cfd_core::compute::gpu::kernels::laplacian::BoundaryCondition) directly."
)]
pub type BoundaryType = BoundaryCondition;
