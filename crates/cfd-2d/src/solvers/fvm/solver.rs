//! FVM solver implementation

use super::config::FvmConfig;
use super::flux::{FluxScheme, FluxSchemeFactory};
use super::geometry::Face;
use crate::fields::Field2D;
use crate::grid::StructuredGrid2D;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// Finite Volume Method solver for 2D problems
pub struct FvmSolver<T: RealField + Copy> {
    config: FvmConfig<T>,
    grid: StructuredGrid2D<T>,
    faces: Vec<Face<T>>,
    flux_scheme: FluxScheme,
}

impl<T: RealField + Copy + FromPrimitive> FvmSolver<T> {
    /// Create a new FVM solver
    pub fn new(config: FvmConfig<T>, grid: StructuredGrid2D<T>) -> Self {
        let faces = Self::generate_faces(&grid);
        Self {
            config,
            grid,
            faces,
            flux_scheme: FluxScheme::Upwind,
        }
    }

    /// Set the flux scheme
    pub fn set_flux_scheme(&mut self, scheme: FluxScheme) {
        self.flux_scheme = scheme;
    }

    /// Generate faces for the structured grid
    fn generate_faces(grid: &StructuredGrid2D<T>) -> Vec<Face<T>> {
        let mut faces = Vec::new();
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;

        // Generate x-direction faces
        for j in 0..ny {
            for i in 0..=nx {
                let center = Vector2::new(
                    T::from_usize(i).unwrap_or_else(T::zero) * dx,
                    (T::from_usize(j).unwrap_or_else(T::zero)
                        + T::from_f64(0.5).unwrap_or_else(T::zero))
                        * dy,
                );
                let normal = Vector2::new(T::one(), T::zero());

                let owner = if i > 0 { Some(j * nx + i - 1) } else { None };
                let neighbor = if i < nx { Some(j * nx + i) } else { None };

                if let Some(o) = owner {
                    faces.push(Face::new(center, normal, dy, o, neighbor));
                }
            }
        }

        // Generate y-direction faces
        for j in 0..=ny {
            for i in 0..nx {
                let center = Vector2::new(
                    (T::from_usize(i).unwrap_or_else(T::zero)
                        + T::from_f64(0.5).unwrap_or_else(T::zero))
                        * dx,
                    T::from_usize(j).unwrap_or_else(T::zero) * dy,
                );
                let normal = Vector2::new(T::zero(), T::one());

                let owner = if j > 0 { Some((j - 1) * nx + i) } else { None };
                let neighbor = if j < ny { Some(j * nx + i) } else { None };

                if let Some(o) = owner {
                    faces.push(Face::new(center, normal, dx, o, neighbor));
                }
            }
        }

        faces
    }

    /// Solve the system
    pub fn solve(
        &mut self,
        phi: &mut Field2D<T>,
        velocity: &Field2D<Vector2<T>>,
        source: &Field2D<T>,
    ) -> Result<()> {
        let flux_calculator = FluxSchemeFactory::create::<T>(self.flux_scheme);

        // Iterative solution
        for _iter in 0..self.config.max_iterations {
            let mut residual = T::zero();

            // Update each cell
            for j in 0..self.config.ny {
                for i in 0..self.config.nx {
                    // Calculate fluxes and update
                    // This is a simplified implementation
                    let phi_old = phi.at(i, j);

                    // Apply relaxation
                    let phi_new = phi_old * (T::one() - self.config.relaxation_factor)
                        + source.at(i, j) * self.config.relaxation_factor;

                    *phi.at_mut(i, j) = phi_new;

                    residual = residual.max((phi_new - phi_old).abs());
                }
            }

            // Check convergence
            if residual < self.config.convergence_tolerance {
                break;
            }
        }

        Ok(())
    }
}
