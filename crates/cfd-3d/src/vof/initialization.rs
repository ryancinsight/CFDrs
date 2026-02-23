//! VOF initialization methods

use super::config::INTERFACE_THICKNESS;
use super::solver::VofSolver;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Initialization methods for VOF
pub enum Initialization<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
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
    },
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> Initialization<T> {
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
            }
            Initialization::Plane { point, normal } => {
                initialize_plane(solver, point, normal);
            }
        }
    }
}

fn initialize_sphere<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    solver: &mut VofSolver<T>,
    center: Vector3<T>,
    radius: T,
) {
    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let x = (T::from_usize(i).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dx;
                let y = (T::from_usize(j).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dy;
                let z = (T::from_usize(k).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dz;

                let pos = Vector3::new(x, y, z);
                let distance = (pos - center).norm();

                let idx = solver.index(i, j, k);

                // Smooth initialization using tanh function
                let eps = <T as FromPrimitive>::from_f64(INTERFACE_THICKNESS).unwrap_or(T::zero()) * solver.dx;
                let arg = (radius - distance) / eps;
                solver.alpha[idx] = <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()) * (T::one() + <T as num_traits::Float>::tanh(arg));
            }
        }
    }
}

fn initialize_block<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    solver: &mut VofSolver<T>,
    min_corner: Vector3<T>,
    max_corner: Vector3<T>,
) {
    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let x = (T::from_usize(i).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dx;
                let y = (T::from_usize(j).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dy;
                let z = (T::from_usize(k).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dz;

                let idx = solver.index(i, j, k);

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
            }
        }
    }
}

fn initialize_plane<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    solver: &mut VofSolver<T>,
    point: Vector3<T>,
    normal: Vector3<T>,
) {
    let normal = normal.normalize();

    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let x = (T::from_usize(i).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dx;
                let y = (T::from_usize(j).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dy;
                let z = (T::from_usize(k).unwrap_or(T::zero())
                    + <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()))
                    * solver.dz;

                let pos = Vector3::new(x, y, z);
                let distance = normal.dot(&(pos - point));

                let idx = solver.index(i, j, k);

                // Smooth initialization
                let eps = <T as FromPrimitive>::from_f64(INTERFACE_THICKNESS).unwrap_or(T::zero()) * solver.dx;
                let arg = distance / eps;
                solver.alpha[idx] = <T as FromPrimitive>::from_f64(0.5).unwrap_or(T::zero()) * (T::one() + <T as num_traits::Float>::tanh(arg));
            }
        }
    }
}
