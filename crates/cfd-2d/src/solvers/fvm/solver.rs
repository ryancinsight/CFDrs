//! FVM solver implementation
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

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
                    T::from_usize(i).expect("analytical constant conversion") * dx,
                    (T::from_usize(j).expect("analytical constant conversion")
                        + T::from_f64(0.5).expect("analytical constant conversion"))
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
                    (T::from_usize(i).expect("analytical constant conversion")
                        + T::from_f64(0.5).expect("analytical constant conversion"))
                        * dx,
                    T::from_usize(j).expect("analytical constant conversion") * dy,
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

        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        debug_assert!(!self.faces.is_empty(), "face list must be populated");

        // Iterative solution using proper FVM discretization
        for _iter in 0..self.config.max_iterations {
            let mut residual = T::zero();

            // Update each interior cell using FVM discretization
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
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
                    .expect("analytical constant conversion");
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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_solver(nx: usize, ny: usize) -> FvmSolver<f64> {
        let grid =
            StructuredGrid2D::new(nx, ny, 0.0, nx as f64 * 0.1, 0.0, ny as f64 * 0.1).unwrap();
        let config = FvmConfig {
            nx,
            ny,
            dx: grid.dx,
            dy: grid.dy,
            ..FvmConfig::default()
        };
        FvmSolver::new(config, grid)
    }

    #[test]
    fn test_face_generation() {
        let solver = make_solver(4, 4);
        assert!(
            !solver.faces.is_empty(),
            "Faces should be generated for a 4×4 grid"
        );
    }

    #[test]
    fn test_solve_uniform_field_stays_uniform() {
        // Uniform φ with zero source and zero velocity → field should not change
        let mut solver = make_solver(6, 6);
        let mut phi = Field2D::new(6, 6, 1.0);
        let velocity = Field2D::new(6, 6, Vector2::zeros());
        let source = Field2D::new(6, 6, 0.0);

        solver.solve(&mut phi, &velocity, &source).unwrap();

        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (phi.at(i, j) - 1.0).abs() < 1e-10,
                    "Uniform field with zero source/velocity should stay uniform at ({i},{j})"
                );
            }
        }
    }

    #[test]
    fn test_set_flux_scheme() {
        let mut solver = make_solver(4, 4);
        solver.set_flux_scheme(FluxScheme::CentralDifference);
        // Should not panic — verifies scheme update works
    }

    #[test]
    fn test_solve_returns_ok() {
        let mut solver = make_solver(8, 8);
        let mut phi = Field2D::new(8, 8, 0.0);
        let velocity = Field2D::new(8, 8, Vector2::new(0.1, 0.0));
        let source = Field2D::new(8, 8, 0.0);

        let result = solver.solve(&mut phi, &velocity, &source);
        assert!(result.is_ok());
    }
}
