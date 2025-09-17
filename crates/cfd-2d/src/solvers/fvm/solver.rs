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
    #[allow(dead_code)]
    grid: StructuredGrid2D<T>,
    #[allow(dead_code)]
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

    /// Solve the system using proper FVM discretization
    pub fn solve(
        &mut self,
        phi: &mut Field2D<T>,
        velocity: &Field2D<Vector2<T>>,
        source: &Field2D<T>,
    ) -> Result<()> {
        let flux_calculator = FluxSchemeFactory::create::<T>(
            self.flux_scheme,
            self.config
                .diffusion_coefficient
                .to_subset()
                .unwrap_or(1e-3),
        );

        let dx = T::from_f64(1.0 / self.config.nx as f64).unwrap_or_else(T::one);
        let dy = T::from_f64(1.0 / self.config.ny as f64).unwrap_or_else(T::one);

        // Iterative solution using proper FVM discretization
        for _iter in 0..self.config.max_iterations {
            let mut residual = T::zero();

            // Update each interior cell using FVM discretization
            for j in 1..self.config.ny - 1 {
                for i in 1..self.config.nx - 1 {
                    let phi_p = phi.at(i, j);
                    let phi_e = phi.at(i + 1, j);
                    let phi_w = phi.at(i - 1, j);
                    let phi_n = phi.at(i, j + 1);
                    let phi_s = phi.at(i, j - 1);

                    // Get velocity components
                    let vel = velocity.at(i, j);
                    let u = vel.x;
                    let v = vel.y;

                    // Calculate convective fluxes using proper discretization
                    let flux_e = flux_calculator.calculate_flux(phi_p, phi_e, phi_p, u, dx)?;
                    let flux_w = flux_calculator.calculate_flux(phi_w, phi_p, phi_w, u, dx)?;
                    let flux_n = flux_calculator.calculate_flux(phi_p, phi_n, phi_p, v, dy)?;
                    let flux_s = flux_calculator.calculate_flux(phi_s, phi_p, phi_s, v, dy)?;

                    // Diffusive terms (central difference)
                    let diff_coeff = T::from_f64(
                        self.config
                            .diffusion_coefficient
                            .to_subset()
                            .unwrap_or(1e-3),
                    )
                    .unwrap_or_else(T::zero);
                    let diff_e = diff_coeff * (phi_e - phi_p) / (dx * dx);
                    let diff_w = diff_coeff * (phi_p - phi_w) / (dx * dx);
                    let diff_n = diff_coeff * (phi_n - phi_p) / (dy * dy);
                    let diff_s = diff_coeff * (phi_p - phi_s) / (dy * dy);

                    // Net flux (convection + diffusion)
                    let net_flux = -(flux_e - flux_w) / dx - (flux_n - flux_s) / dy
                        + (diff_e - diff_w)
                        + (diff_n - diff_s);

                    // Update with source term and relaxation
                    let phi_new =
                        phi_p + self.config.relaxation_factor * (net_flux + source.at(i, j));

                    if let Some(p) = phi.at_mut(i, j) {
                        *p = phi_new;
                    }

                    residual = residual.max((phi_new - phi_p).abs());
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
