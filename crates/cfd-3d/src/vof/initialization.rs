//! VOF initialization methods
//!
//! # Theorem — Signed-Distance VOF Initialisation
//!
//! When initialising the volume fraction field from a signed-distance
//! function $\phi(\mathbf{x})$, the smoothed Heaviside
//!
//! ```text
//! α(x) = H_ε(φ(x)) = ½(1 + φ/ε + sin(πφ/ε)/π)  for |φ| ≤ ε
//! ```
//!
//! with interface thickness $\varepsilon = O(h)$ preserves the enclosed volume
//! to $O(h^2)$ on a uniform grid of spacing $h$.

use super::config::INTERFACE_THICKNESS;
use super::scalar::{self, VofScalar};
use super::solver::VofSolver;
use leto::geometry::Vector3;

/// Initialization methods for VOF
pub enum Initialization<T: VofScalar> {
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

impl<T: VofScalar> Initialization<T> {
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

fn initialize_sphere<T: VofScalar>(solver: &mut VofSolver<T>, center: Vector3<T>, radius: T) {
    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let half = scalar::constant::<T>(0.5);
                let x = (scalar::from_usize::<T>(i) + half) * solver.dx;
                let y = (scalar::from_usize::<T>(j) + half) * solver.dy;
                let z = (scalar::from_usize::<T>(k) + half) * solver.dz;

                let pos = Vector3::new(x, y, z);
                let distance = (pos - center).norm();

                let idx = solver.index(i, j, k);

                // Smooth initialization using tanh function
                let eps = scalar::constant::<T>(INTERFACE_THICKNESS) * solver.dx;
                let arg = (radius - distance) / eps;
                solver.alpha[idx] = half * (scalar::one::<T>() + scalar::tanh(arg));
            }
        }
    }
}

fn initialize_block<T: VofScalar>(
    solver: &mut VofSolver<T>,
    min_corner: Vector3<T>,
    max_corner: Vector3<T>,
) {
    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let half = scalar::constant::<T>(0.5);
                let x = (scalar::from_usize::<T>(i) + half) * solver.dx;
                let y = (scalar::from_usize::<T>(j) + half) * solver.dy;
                let z = (scalar::from_usize::<T>(k) + half) * solver.dz;

                let idx = solver.index(i, j, k);

                // Check if point is inside block
                if x >= min_corner.x
                    && x <= max_corner.x
                    && y >= min_corner.y
                    && y <= max_corner.y
                    && z >= min_corner.z
                    && z <= max_corner.z
                {
                    solver.alpha[idx] = scalar::one();
                } else {
                    solver.alpha[idx] = scalar::zero();
                }
            }
        }
    }
}

fn initialize_plane<T: VofScalar>(
    solver: &mut VofSolver<T>,
    point: Vector3<T>,
    normal: Vector3<T>,
) {
    let normal = normal.normalize();

    for k in 0..solver.nz {
        for j in 0..solver.ny {
            for i in 0..solver.nx {
                let half = scalar::constant::<T>(0.5);
                let x = (scalar::from_usize::<T>(i) + half) * solver.dx;
                let y = (scalar::from_usize::<T>(j) + half) * solver.dy;
                let z = (scalar::from_usize::<T>(k) + half) * solver.dz;

                let pos = Vector3::new(x, y, z);
                let distance = normal.dot(pos - point);

                let idx = solver.index(i, j, k);

                // Smooth initialization
                let eps = scalar::constant::<T>(INTERFACE_THICKNESS) * solver.dx;
                let arg = distance / eps;
                solver.alpha[idx] = half * (scalar::one::<T>() + scalar::tanh(arg));
            }
        }
    }
}
