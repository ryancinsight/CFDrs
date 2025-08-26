//! Flow field operations
//!
//! Provides operations on flow fields including vorticity, divergence,
//! and other fluid mechanical quantities.

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use rayon::prelude::*;
use super::fields::VelocityField;
/// Operations on flow fields
pub struct FlowOperations;
impl FlowOperations {
    /// Calculate vorticity field (curl of velocity)
    pub fn vorticity<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity: &VelocityField<T>,
    ) -> Vec<Vector3<T>> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut vorticity = vec![Vector3::zeros(); nx * ny * nz];
        // Use parallel iteration for better performance
        vorticity
            .par_iter_mut()
            .enumerate()
            .for_each(|(idx, vort)| {
                let k = idx / (nx * ny);
                let j = (idx % (nx * ny)) / nx;
                let i = idx % nx;
                *vort = vorticity_at_point(velocity, i, j, k, nx, ny, nz);
            });
    }
    /// Calculate divergence field
    pub fn divergence<T: RealField + Copy + FromPrimitive + Send + Sync>(
    ) -> Vec<T> {
        let mut divergence = vec![T::zero(); nx * ny * nz];
        divergence
            .for_each(|(idx, div)| {
                let dx = T::one();
                let dy = T::one();
                let dz = T::one();
                // Central differences with boundary handling
                let dudx = if i > 0 && i < nx - 1 {
                    let idx_plus = k * nx * ny + j * nx + (i + 1);
                    let idx_minus = k * nx * ny + j * nx + (i - 1);
                    (velocity.components[idx_plus].x - velocity.components[idx_minus].x)
                        / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dx)
                } else {
                    T::zero()
                };
                let dvdy = if j > 0 && j < ny - 1 {
                    let idx_plus = k * nx * ny + (j + 1) * nx + i;
                    let idx_minus = k * nx * ny + (j - 1) * nx + i;
                    (velocity.components[idx_plus].y - velocity.components[idx_minus].y)
                        / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dy)
                let dwdz = if k > 0 && k < nz - 1 {
                    let idx_plus = (k + 1) * nx * ny + j * nx + i;
                    let idx_minus = (k - 1) * nx * ny + j * nx + i;
                    (velocity.components[idx_plus].z - velocity.components[idx_minus].z)
                        / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dz)
                *div = dudx + dvdy + dwdz;
    /// Calculate kinetic energy field
    pub fn kinetic_energy<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity
            .components
            .par_iter()
            .map(|v| {
                T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()))
                    * v.norm_squared()
            })
            .collect()
    /// Calculate enstrophy field
    pub fn enstrophy<T: RealField + Copy + FromPrimitive + Send + Sync>(
        let vorticity = Self::vorticity(velocity);
            .map(|w| {
                    * w.norm_squared()
}
// Helper function for vorticity calculation at a point
fn vorticity_at_point<T: RealField + Copy + FromPrimitive>(
    velocity: &VelocityField<T>,
    i: usize,
    j: usize,
    k: usize,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Vector3<T> {
    let dx = T::one();
    let dy = T::one();
    let dz = T::one();
    // Calculate velocity gradients using central differences
    let dvdz = if k > 0 && k < nz - 1 {
        let idx_plus = (k + 1) * nx * ny + j * nx + i;
        let idx_minus = (k - 1) * nx * ny + j * nx + i;
        (velocity.components[idx_plus].y - velocity.components[idx_minus].y)
            / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dz)
    } else {
        T::zero()
    };
    let dwdy = if j > 0 && j < ny - 1 {
        let idx_plus = k * nx * ny + (j + 1) * nx + i;
        let idx_minus = k * nx * ny + (j - 1) * nx + i;
        (velocity.components[idx_plus].z - velocity.components[idx_minus].z)
            / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dy)
    let dudz = if k > 0 && k < nz - 1 {
        (velocity.components[idx_plus].x - velocity.components[idx_minus].x)
    let dwdx = if i > 0 && i < nx - 1 {
        let idx_plus = k * nx * ny + j * nx + (i + 1);
        let idx_minus = k * nx * ny + j * nx + (i - 1);
            / (T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()) * dx)
    let dvdx = if i > 0 && i < nx - 1 {
    let dudy = if j > 0 && j < ny - 1 {
    // Vorticity = curl(velocity)
    Vector3::new(dwdy - dvdz, dudz - dwdx, dvdx - dudy)
