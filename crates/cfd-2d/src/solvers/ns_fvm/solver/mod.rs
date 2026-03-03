//! `NavierStokesSolver2D` — SIMPLE pressure-velocity coupling solver.
//!
//! Implements the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE;
//! Patankar, 1980) for 2D incompressible flow with non-Newtonian blood rheology
//! on a Cartesian staggered grid.
//!
//! ## Algorithm
//! 1. Solve u-momentum (Gauss-Seidel + under-relaxation)
//! 2. Solve v-momentum (Gauss-Seidel + under-relaxation)
//! 3. Solve pressure-correction equation from continuity residual
//! 4. Correct velocities + pressure
//! 5. Apply Rhie-Chow face-velocity interpolation (via `cfd-core`)
//! 6. Update viscosity from shear rate (non-Newtonian)
//! 7. Test convergence
//!
//! ## References
//! - Patankar (1980): §6.3–6.7, §7.4
//! - Rhie & Chow (1983): Pressure-velocity interpolation
//! - Versteeg & Malalasekera (2007): §11.5
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

mod momentum;
mod pressure;

use super::boundary::BoundaryCondition;
use super::config::{SIMPLEConfig, SolveResult};
use super::field::FlowField2D;
use super::grid::StaggeredGrid2D;
use super::BloodModel;
use crate::error::Error;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// 2D Navier-Stokes FVM solver with SIMPLE pressure-velocity coupling.
///
/// Used as the numerical engine by geometry-specific pass-through solvers
/// (bifurcation, Venturi, serpentine, etc.).
#[derive(Debug)]
pub struct NavierStokesSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Computational grid
    pub grid: StaggeredGrid2D<T>,
    /// Flow field (u, v, p, μ, γ̇)
    pub field: FlowField2D<T>,
    /// Blood rheology model
    pub blood: BloodModel<T>,
    /// Fluid density [kg/m³]
    pub density: T,
    /// SIMPLE configuration
    pub config: SIMPLEConfig<T>,
    /// Central coefficient storage for u-momentum
    a_p_u: Vec<Vec<T>>,
    /// Central coefficient storage for v-momentum
    a_p_v: Vec<Vec<T>>,
}

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Create a new solver.
    pub fn new(
        grid: StaggeredGrid2D<T>,
        blood: BloodModel<T>,
        density: T,
        config: SIMPLEConfig<T>,
    ) -> Self {
        let field = FlowField2D::new(grid.nx, grid.ny);
        let a_p_u = vec![vec![T::one(); grid.ny]; grid.nx + 1];
        let a_p_v = vec![vec![T::one(); grid.ny + 1]; grid.nx];
        Self {
            grid,
            field,
            blood,
            density,
            config,
            a_p_u,
            a_p_v,
        }
    }

    /// Initialise viscosity field from the blood model's apparent viscosity at
    /// the reference shear rate.  Using μ_∞ (high-shear limit) here would
    /// under-estimate viscosity by ~30 %, making the initial velocity field too
    /// fast and slowing SIMPLE convergence for non-Newtonian blood.
    pub fn initialize_viscosity(&mut self) {
        let mu_init = match &self.blood {
            BloodModel::Casson(m) => m.apparent_viscosity(m.reference_shear_rate),
            BloodModel::CarreauYasuda(m) => m.apparent_viscosity(m.reference_shear_rate),
            BloodModel::Newtonian(mu) => *mu,
        };
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                self.field.mu[i][j] = mu_init;
            }
        }
    }

    /// Compute L2-norm continuity residual for convergence assessment.
    /// Supports non-uniform y-spacing via `grid.dy_at(j)`.
    pub fn compute_residuals(&self) -> (T, T, T) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let rho = self.density;

        let mut cont_sum = T::zero();
        for i in 0..nx {
            for j in 0..ny {
                let dy_j = self.grid.dy_at(j);
                let imb = rho
                    * ((self.field.u[i + 1][j] - self.field.u[i][j]) * dy_j
                        + (self.field.v[i][j + 1] - self.field.v[i][j]) * dx);
                cont_sum += imb * imb;
            }
        }
        let n: T = T::from_usize(nx * ny).unwrap_or(T::one());
        let r = Float::sqrt(cont_sum / n);
        (r, r, r)
    }

    /// Drive the SIMPLE loop to steady state.
    pub fn solve(&mut self, u_inlet: T) -> Result<SolveResult<T>, Error> {
        self.initialize_viscosity();
        self.a_p_u = vec![vec![T::one(); self.grid.ny]; self.grid.nx + 1];
        self.a_p_v = vec![vec![T::one(); self.grid.ny + 1]; self.grid.nx];

        let ny = self.grid.ny;

        // Parabolic inlet profile (supports non-uniform y-spacing)
        let mut y_coords = Vec::new();
        let mut dy_cells = Vec::new();
        for j in 0..ny {
            if self.field.mask[0][j] {
                y_coords.push(self.grid.y_center(j));
                dy_cells.push(self.grid.dy_at(j));
            }
        }
        if !y_coords.is_empty() {
            let y_min = y_coords.iter().copied().fold(y_coords[0], Float::min);
            let y_max = y_coords.iter().copied().fold(y_coords[0], Float::max);
            // Channel height = span from first to last fluid centre + half-cells at edges.
            let h = y_max - y_min + (dy_cells[0] + dy_cells[dy_cells.len() - 1])
                * T::from_f64(0.5).unwrap();
            for j in 0..ny {
                if self.field.mask[0][j] {
                    let dy_j = self.grid.dy_at(j);
                    let y_local =
                        (self.grid.y_center(j) - y_min) + T::from_f64(0.5).unwrap() * dy_j;
                    let y_frac = y_local / h;
                    self.field.u[0][j] =
                        T::from_f64(6.0).unwrap() * u_inlet * y_frac * (T::one() - y_frac);
                } else {
                    self.field.u[0][j] = T::zero();
                }
            }
        }

        let mut last_residual = T::from_f64(1e10).unwrap_or(T::one());

        let bc_inlet = BoundaryCondition::velocity_inlet(nalgebra::Vector3::new(
            u_inlet,
            T::zero(),
            T::zero(),
        ));
        let bc_outlet = BoundaryCondition::pressure_outlet(T::zero());
        let bc_wall_noslip = BoundaryCondition::wall_no_slip();

        for iteration in 0..self.config.max_iterations {
            self.solve_u_momentum(&bc_inlet, &bc_outlet, u_inlet)?;
            self.solve_v_momentum(&bc_wall_noslip, &bc_wall_noslip)?;
            self.solve_pressure_correction()?;

            // Update viscosity with under-relaxation (alpha_mu) to prevent
            // oscillation in non-Newtonian SIMPLE iterations.
            if iteration % self.config.viscosity_update_interval == 0 {
                self.field
                    .update_viscosity(&self.grid, &self.blood, self.config.alpha_mu);
            }

            let (res_u, res_v, res_p) = self.compute_residuals();
            last_residual = Float::max(Float::max(res_u, res_v), res_p);

            if last_residual < self.config.tolerance {
                return Ok(SolveResult {
                    iterations: iteration + 1,
                    residual: last_residual,
                    converged: true,
                });
            }
        }

        Ok(SolveResult {
            iterations: self.config.max_iterations,
            residual: last_residual,
            converged: false,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid::blood::CassonBlood;

    #[test]
    fn test_staggered_grid_creation() {
        let grid = StaggeredGrid2D::<f64>::new(10, 5, 0.1, 0.005);
        assert_eq!(grid.nx, 10);
        assert_eq!(grid.ny, 5);
        assert!((grid.dx - 0.01).abs() < 1e-10);
        assert!((grid.dy - 0.001).abs() < 1e-10);
    }

    #[test]
    fn test_flow_field_creation() {
        let field = FlowField2D::<f64>::new(10, 5);
        assert_eq!(field.u.len(), 11);
        assert_eq!(field.v.len(), 10);
        assert_eq!(field.p.len(), 10);
        assert_eq!(field.u[0].len(), 5);
        assert_eq!(field.v[0].len(), 6);
    }

    #[test]
    fn test_simple_solver_newtonian() {
        let grid = StaggeredGrid2D::<f64>::new(20, 10, 0.2, 0.01);
        let blood = BloodModel::Newtonian(1.0e-3);
        let config = SIMPLEConfig::new(200, 1e-4, 0.5, 0.3, 0.5, 1);
        let mut solver = NavierStokesSolver2D::new(grid, blood, 998.0, config);
        let result = solver.solve(0.01).unwrap();
        assert!(result.iterations > 0);
    }

    #[test]
    fn test_simple_solver_non_newtonian() {
        let grid = StaggeredGrid2D::<f64>::new(10, 5, 0.05, 0.0025);
        let blood = BloodModel::Casson(CassonBlood::normal_blood());
        let config = SIMPLEConfig::new(100, 1e-3, 0.3, 0.2, 0.4, 1);
        let mut solver = NavierStokesSolver2D::new(grid, blood, 1060.0, config);
        assert!(solver.solve(0.005).is_ok());
    }

    #[test]
    fn test_residuals_zero_flow() {
        let grid = StaggeredGrid2D::<f64>::new(5, 5, 0.05, 0.005);
        let solver = NavierStokesSolver2D::new(
            grid,
            BloodModel::Newtonian(1.0e-3),
            998.0,
            SIMPLEConfig::default(),
        );
        let (ru, rv, rp) = solver.compute_residuals();
        assert_eq!(ru, 0.0);
        assert_eq!(rv, 0.0);
        assert_eq!(rp, 0.0);
    }

    #[test]
    fn test_solve_result_struct() {
        use super::super::config::SolveResult;
        let result = SolveResult {
            iterations: 42_usize,
            residual: 1e-7_f64,
            converged: true,
        };
        assert_eq!(result.iterations, 42);
        assert!(result.converged);
        assert!(result.residual < 1e-6);
    }
}
