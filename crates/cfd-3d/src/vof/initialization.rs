//! VOF initialization methods

use super::config::INTERFACE_THICKNESS;
use cfd_core::numeric;
use super::solver::VofSolver;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
/// Initialization methods for VOF
pub enum Initialization<T: RealField + Copy> {
    /// Initialize with a sphere
    Sphere { center: Vector3<T>, radius: T },
    /// Initialize with a rectangular block
    Block {
        min_corner: Vector3<T>,
        max_corner: Vector3<T>,
    },
    /// Initialize with a plane
    Plane {
        point: Vector3<T>,
        normal: Vector3<T>,
}
impl<T: RealField + FromPrimitive + Copy> Initialization<T> {
    /// Apply initialization to the solver
    pub fn apply(self, solver: &mut VofSolver<T>) {
        match self {
            Initialization::Sphere { center, radius } => {
                initialize_sphere(solver, center, radius);
            }
            Initialization::Block {
                min_corner,
                max_corner,
            } => {
                initialize_block(solver, min_corner, max_corner);
            Initialization::Plane { point, normal } => {
                initialize_plane(solver, point, normal);
        }
    }
fn initialize_sphere<T: RealField + FromPrimitive + Copy>(
    solver: &mut VofSolver<T>,
    center: Vector3<T>,
    radius: T,
) {
    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let x = (cfd_core::numeric::from_usize(i)?
                    + cfd_core::numeric::from_f64(0.5)?)
                    * solver.dx;
                let y = (cfd_core::numeric::from_usize(j)?
                    * solver.dy;
                let z = (cfd_core::numeric::from_usize(k)?
                    * solver.dz;
                let pos = Vector3::new(x, y, z);
                let distance = (pos - center).norm();
                let idx = solver.index(i, j, k);
                // Smooth initialization using tanh function
                let eps = cfd_core::numeric::from_f64(INTERFACE_THICKNESS)? * solver.dx;
                let arg = (radius - distance) / eps;
                solver.alpha[idx] = cfd_core::numeric::from_f64(0.5)? * (T::one() + arg.tanh());
fn initialize_block<T: RealField + FromPrimitive + Copy>(
    min_corner: Vector3<T>,
    max_corner: Vector3<T>,
                // Check if point is inside block
                if x >= min_corner.x
                    && x <= max_corner.x
                    && y >= min_corner.y
                    && y <= max_corner.y
                    && z >= min_corner.z
                    && z <= max_corner.z
                {
                    solver.alpha[idx] = T::one();
                } else {
                    solver.alpha[idx] = T::zero();
                }
fn initialize_plane<T: RealField + FromPrimitive + Copy>(
    point: Vector3<T>,
    normal: Vector3<T>,
    let normal = normal.normalize();
                let distance = normal.dot(&(pos - point));
                // Smooth initialization
                let arg = distance / eps;
