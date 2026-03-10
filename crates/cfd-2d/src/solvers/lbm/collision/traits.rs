//! Collision operator trait for LBM solvers.
//!
//! All collision operators work on the flat distribution buffer
//! with layout `f[j * nx * 9 + i * 9 + q]`.

use nalgebra::RealField;

/// Trait for LBM collision operators operating on flat distribution buffers.
pub trait CollisionOperator<T: RealField + Copy>: Send + Sync {
    /// Apply collision to the flat distribution buffer `f`.
    ///
    /// # Arguments
    /// * `f`        - flat distribution slice, layout `f[j*nx*9 + i*9 + q]`
    /// * `density`  - flat density slice, layout `density[j*nx + i]`
    /// * `velocity` - flat velocity slice, layout `velocity[(j*nx + i)*2 + d]`
    /// * `nx`, `ny` - grid dimensions
    fn collide(&self, f: &mut [T], density: &[T], velocity: &[T], nx: usize, ny: usize);

    /// Return the primary relaxation time τ.
    fn tau(&self) -> T;

    /// Return the kinematic viscosity ν = c_s²(τ − ½)Δt / Δx² · Δx².
    fn viscosity(&self, dt: T, dx: T) -> T;
}
