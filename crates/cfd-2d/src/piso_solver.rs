//! Refactored PISO solver with clean separation of concerns
//!
//! This implementation follows SOLID principles with proper separation of:
//! - Configuration (immutable solver parameters)
//! - State (persistent simulation data)
//! - Workspace (temporary per-timestep data)
//! - Discretization schemes (configurable numerical methods)

use crate::field::{ScalarField, VectorField, SolverState, PisoWorkspace};
use crate::grid::StructuredGrid2D;
use cfd_core::{Result, BoundaryCondition, constants};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Discretization scheme for gradient calculations
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum GradientScheme {
    /// Second-order central difference
    CentralDifference,
    /// First-order forward difference
    ForwardDifference,
    /// First-order backward difference
    BackwardDifference,
}

/// Discretization scheme for convection terms
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ConvectionScheme {
    /// First-order upwind
    Upwind,
    /// Hybrid upwind/central
    Hybrid,
    /// Quadratic Upstream Interpolation for Convective Kinematics
    QUICK,
    /// Central difference (can be unstable at high Pe)
    Central,
}

/// PISO solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PisoConfig<T: RealField> {
    /// Time step for transient simulation
    pub time_step: T,
    /// Number of corrector steps (typically 2)
    pub num_correctors: usize,
    /// Convergence tolerance for velocity
    pub velocity_tolerance: T,
    /// Convergence tolerance for pressure
    pub pressure_tolerance: T,
    /// Number of non-orthogonal correctors
    pub non_orthogonal_correctors: usize,
    /// Maximum iterations for pressure correction
    pub pressure_correction_max_iters: usize,
    /// Maximum iterations for momentum predictor
    pub momentum_max_iters: usize,
    /// Gradient discretization scheme
    pub gradient_scheme: GradientScheme,
    /// Convection discretization scheme
    pub convection_scheme: ConvectionScheme,
}

impl<T: RealField + FromPrimitive> PisoConfig<T> {
    /// Create a new configuration with validated parameters
    pub fn new() -> Result<Self> {
        Ok(Self {
            time_step: T::from_f64(0.01)
                .ok_or_else(|| cfd_core::Error::InvalidConfiguration(
                    "Cannot convert time_step to target type".into()
                ))?,
            num_correctors: 2,
            velocity_tolerance: T::from_f64(1e-6)
                .ok_or_else(|| cfd_core::Error::InvalidConfiguration(
                    "Cannot convert velocity_tolerance to target type".into()
                ))?,
            pressure_tolerance: T::from_f64(1e-6)
                .ok_or_else(|| cfd_core::Error::InvalidConfiguration(
                    "Cannot convert pressure_tolerance to target type".into()
                ))?,
            non_orthogonal_correctors: 1,
            pressure_correction_max_iters: 20,
            momentum_max_iters: 50,
            gradient_scheme: GradientScheme::CentralDifference,
            convection_scheme: ConvectionScheme::Hybrid,
        })
    }

    /// Builder pattern for configuration
    pub fn builder() -> PisoConfigBuilder<T> {
        PisoConfigBuilder::default()
    }
}

/// Builder for PISO configuration
#[derive(Debug, Clone)]
pub struct PisoConfigBuilder<T: RealField> {
    time_step: Option<T>,
    num_correctors: Option<usize>,
    velocity_tolerance: Option<T>,
    pressure_tolerance: Option<T>,
    non_orthogonal_correctors: Option<usize>,
    pressure_correction_max_iters: Option<usize>,
    momentum_max_iters: Option<usize>,
    gradient_scheme: Option<GradientScheme>,
    convection_scheme: Option<ConvectionScheme>,
}

impl<T: RealField> Default for PisoConfigBuilder<T> {
    fn default() -> Self {
        Self {
            time_step: None,
            num_correctors: None,
            velocity_tolerance: None,
            pressure_tolerance: None,
            non_orthogonal_correctors: None,
            pressure_correction_max_iters: None,
            momentum_max_iters: None,
            gradient_scheme: None,
            convection_scheme: None,
        }
    }
}

impl<T: RealField + FromPrimitive> PisoConfigBuilder<T> {
    pub fn time_step(mut self, dt: T) -> Self {
        self.time_step = Some(dt);
        self
    }

    pub fn num_correctors(mut self, n: usize) -> Self {
        self.num_correctors = Some(n);
        self
    }

    pub fn velocity_tolerance(mut self, tol: T) -> Self {
        self.velocity_tolerance = Some(tol);
        self
    }

    pub fn pressure_tolerance(mut self, tol: T) -> Self {
        self.pressure_tolerance = Some(tol);
        self
    }

    pub fn gradient_scheme(mut self, scheme: GradientScheme) -> Self {
        self.gradient_scheme = Some(scheme);
        self
    }

    pub fn convection_scheme(mut self, scheme: ConvectionScheme) -> Self {
        self.convection_scheme = Some(scheme);
        self
    }

    pub fn build(self) -> Result<PisoConfig<T>> {
        let default = PisoConfig::new()?;
        Ok(PisoConfig {
            time_step: self.time_step.unwrap_or(default.time_step),
            num_correctors: self.num_correctors.unwrap_or(default.num_correctors),
            velocity_tolerance: self.velocity_tolerance.unwrap_or(default.velocity_tolerance),
            pressure_tolerance: self.pressure_tolerance.unwrap_or(default.pressure_tolerance),
            non_orthogonal_correctors: self.non_orthogonal_correctors
                .unwrap_or(default.non_orthogonal_correctors),
            pressure_correction_max_iters: self.pressure_correction_max_iters
                .unwrap_or(default.pressure_correction_max_iters),
            momentum_max_iters: self.momentum_max_iters
                .unwrap_or(default.momentum_max_iters),
            gradient_scheme: self.gradient_scheme.unwrap_or(default.gradient_scheme),
            convection_scheme: self.convection_scheme.unwrap_or(default.convection_scheme),
        })
    }
}

/// PISO solver with clean separation of concerns
pub struct PisoSolver<T: RealField> {
    /// Solver configuration (immutable)
    config: PisoConfig<T>,
    /// Grid structure
    grid: StructuredGrid2D<T>,
    /// Fluid properties
    density: T,
    viscosity: T,
    /// Boundary conditions
    boundary_conditions: HashMap<(usize, usize), BoundaryCondition<T>>,
}

impl<T: RealField + FromPrimitive + Copy> PisoSolver<T> {
    /// Create a new PISO solver
    pub fn new(
        config: PisoConfig<T>,
        grid: StructuredGrid2D<T>,
        density: T,
        viscosity: T,
    ) -> Result<Self> {
        Ok(Self {
            config,
            grid,
            density,
            viscosity,
            boundary_conditions: HashMap::new(),
        })
    }

    /// Set boundary conditions
    pub fn set_boundary_conditions(&mut self, bc: HashMap<(usize, usize), BoundaryCondition<T>>) {
        self.boundary_conditions = bc;
    }

    /// Execute one PISO time step
    pub fn step(&self, state: &mut SolverState<T>) -> Result<()> {
        // Create workspace for this time step
        let mut workspace = PisoWorkspace::new(self.grid.nx, self.grid.ny);
        
        // Swap old and current velocity fields
        state.swap_velocities();
        
        // Step 1: Momentum predictor
        self.solve_momentum_predictor(state, &mut workspace)?;
        
        // Step 2: Pressure corrections
        for _ in 0..self.config.num_correctors {
            self.solve_pressure_correction(state, &mut workspace)?;
            self.correct_velocity(state, &workspace)?;
        }
        
        // Clear workspace (it will be deallocated when it goes out of scope)
        workspace.clear();
        
        Ok(())
    }

    /// Solve momentum predictor step
    fn solve_momentum_predictor(
        &self,
        state: &mut SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        // Compute momentum coefficients
        self.compute_momentum_coefficients(state, workspace)?;
        
        // Build source terms
        self.build_momentum_sources(state, workspace)?;
        
        // Solve momentum equations using Gauss-Seidel
        for _ in 0..self.config.momentum_max_iters {
            self.gauss_seidel_momentum(state, workspace)?;
        }
        
        Ok(())
    }

    /// Compute momentum equation coefficients
    fn compute_momentum_coefficients(
        &self,
        state: &SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        let dt = self.config.time_step;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        
        // Diffusion coefficients
        let diff_x = self.viscosity / dx;
        let diff_y = self.viscosity / dy;
        
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Face velocities (simple averaging)
                let u_e = (state.velocity[(i, j)].x + state.velocity[(i + 1, j)].x) / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                let u_w = (state.velocity[(i - 1, j)].x + state.velocity[(i, j)].x) / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                let v_n = (state.velocity[(i, j)].y + state.velocity[(i, j + 1)].y) / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                let v_s = (state.velocity[(i, j - 1)].y + state.velocity[(i, j)].y) / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                
                // Convective fluxes
                let conv_e = self.density * u_e / dx;
                let conv_w = self.density * u_w / dx;
                let conv_n = self.density * v_n / dy;
                let conv_s = self.density * v_s / dy;
                
                // Apply convection scheme
                match self.config.convection_scheme {
                    ConvectionScheme::Upwind => {
                        workspace.a_e[(i, j)] = diff_x + T::max(T::zero(), -conv_e);
                        workspace.a_w[(i, j)] = diff_x + T::max(T::zero(), conv_w);
                        workspace.a_n[(i, j)] = diff_y + T::max(T::zero(), -conv_n);
                        workspace.a_s[(i, j)] = diff_y + T::max(T::zero(), conv_s);
                    }
                    ConvectionScheme::Hybrid => {
                        let pe_x = conv_e / diff_x;
                        let pe_y = conv_n / diff_y;
                        
                        if pe_x.abs() < T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? {
                            // Central difference
                            workspace.a_e[(i, j)] = diff_x - conv_e / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                            workspace.a_w[(i, j)] = diff_x + conv_w / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                        } else {
                            // Upwind
                            workspace.a_e[(i, j)] = diff_x + T::max(T::zero(), -conv_e);
                            workspace.a_w[(i, j)] = diff_x + T::max(T::zero(), conv_w);
                        }
                        
                        if pe_y.abs() < T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? {
                            workspace.a_n[(i, j)] = diff_y - conv_n / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                            workspace.a_s[(i, j)] = diff_y + conv_s / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                        } else {
                            workspace.a_n[(i, j)] = diff_y + T::max(T::zero(), -conv_n);
                            workspace.a_s[(i, j)] = diff_y + T::max(T::zero(), conv_s);
                        }
                    }
                    ConvectionScheme::Central => {
                        workspace.a_e[(i, j)] = diff_x - conv_e / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                        workspace.a_w[(i, j)] = diff_x + conv_w / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                        workspace.a_n[(i, j)] = diff_y - conv_n / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                        workspace.a_s[(i, j)] = diff_y + conv_s / T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                    }
                    ConvectionScheme::QUICK => {
                        // Simplified QUICK implementation
                        // Full QUICK would require more stencil points
                        workspace.a_e[(i, j)] = diff_x + T::max(T::zero(), -conv_e);
                        workspace.a_w[(i, j)] = diff_x + T::max(T::zero(), conv_w);
                        workspace.a_n[(i, j)] = diff_y + T::max(T::zero(), -conv_n);
                        workspace.a_s[(i, j)] = diff_y + T::max(T::zero(), conv_s);
                    }
                }
                
                // Central coefficient (includes transient term)
                workspace.a_p[(i, j)] = workspace.a_e[(i, j)] + workspace.a_w[(i, j)]
                    + workspace.a_n[(i, j)] + workspace.a_s[(i, j)]
                    + self.density * dx * dy / dt;
            }
        }
        
        Ok(())
    }

    /// Build momentum equation source terms
    fn build_momentum_sources(
        &self,
        state: &SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        let dt = self.config.time_step;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Pressure gradient
                let (grad_p_x, grad_p_y) = self.compute_pressure_gradient(state, i, j)?;
                
                // Transient term
                let cell_volume = dx * dy;
                workspace.b_x[(i, j)] = self.density * cell_volume * state.velocity_old[(i, j)].x / dt
                    - grad_p_x * cell_volume;
                workspace.b_y[(i, j)] = self.density * cell_volume * state.velocity_old[(i, j)].y / dt
                    - grad_p_y * cell_volume;
            }
        }
        
        Ok(())
    }

    /// Compute pressure gradient at a cell
    fn compute_pressure_gradient(
        &self,
        state: &SolverState<T>,
        i: usize,
        j: usize,
    ) -> Result<(T, T)> {
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        
        let (grad_x, grad_y) = match self.config.gradient_scheme {
            GradientScheme::CentralDifference => {
                let two = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
                (
                    (state.pressure[(i + 1, j)] - state.pressure[(i - 1, j)]) / (two * dx),
                    (state.pressure[(i, j + 1)] - state.pressure[(i, j - 1)]) / (two * dy),
                )
            }
            GradientScheme::ForwardDifference => {
                (
                    (state.pressure[(i + 1, j)] - state.pressure[(i, j)]) / dx,
                    (state.pressure[(i, j + 1)] - state.pressure[(i, j)]) / dy,
                )
            }
            GradientScheme::BackwardDifference => {
                (
                    (state.pressure[(i, j)] - state.pressure[(i - 1, j)]) / dx,
                    (state.pressure[(i, j)] - state.pressure[(i, j - 1)]) / dy,
                )
            }
        };
        
        Ok((grad_x, grad_y))
    }

    /// Gauss-Seidel iteration for momentum equations
    fn gauss_seidel_momentum(
        &self,
        state: &mut SolverState<T>,
        workspace: &PisoWorkspace<T>,
    ) -> Result<()> {
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Skip if diagonal coefficient is too small
                if workspace.a_p[(i, j)].abs() <= T::from_f64(1e-14).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? {
                    continue;
                }
                
                // X-momentum equation
                let u_x_neighbors = workspace.a_w[(i, j)] * state.velocity[(i - 1, j)].x
                    + workspace.a_e[(i, j)] * state.velocity[(i + 1, j)].x
                    + workspace.a_s[(i, j)] * state.velocity[(i, j - 1)].x
                    + workspace.a_n[(i, j)] * state.velocity[(i, j + 1)].x;
                
                state.velocity[(i, j)].x = (workspace.b_x[(i, j)] - u_x_neighbors) / workspace.a_p[(i, j)];
                
                // Y-momentum equation
                let u_y_neighbors = workspace.a_w[(i, j)] * state.velocity[(i - 1, j)].y
                    + workspace.a_e[(i, j)] * state.velocity[(i + 1, j)].y
                    + workspace.a_s[(i, j)] * state.velocity[(i, j - 1)].y
                    + workspace.a_n[(i, j)] * state.velocity[(i, j + 1)].y;
                
                state.velocity[(i, j)].y = (workspace.b_y[(i, j)] - u_y_neighbors) / workspace.a_p[(i, j)];
            }
        }
        
        // Apply boundary conditions
        self.apply_velocity_boundary_conditions(state)?;
        
        Ok(())
    }

    /// Solve pressure correction equation
    /// Solves: ∇²p' = (ρ/Δt)∇·u* using Jacobi iteration
    /// Reference: Issa (1986), J. Comp. Phys. 62, 40-65
    fn solve_pressure_correction(
        &self,
        state: &SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let dt = self.config.time_step;
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let factor = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * (dx2 + dy2);
        
        // Build RHS: (ρ/Δt)∇·u*
        let mut rhs = ScalarField::new(self.grid.nx, self.grid.ny);
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let div_u = self.compute_velocity_divergence(state, i, j)?;
                rhs[(i, j)] = self.density * div_u / dt;
            }
        }
        
        // Solve Poisson equation using Jacobi iteration
        for _ in 0..self.config.pressure_correction_max_iters {
            let p_old = workspace.p_prime.clone();
            
            for i in 1..self.grid.nx - 1 {
                for j in 1..self.grid.ny - 1 {
                    // Laplacian stencil: (p_E + p_W)/dx² + (p_N + p_S)/dy² - rhs
                    let laplacian_x = (p_old[(i+1, j)] + p_old[(i-1, j)]) * dy2;
                    let laplacian_y = (p_old[(i, j+1)] + p_old[(i, j-1)]) * dx2;
                    workspace.p_prime[(i, j)] = (laplacian_x + laplacian_y - rhs[(i, j)] * dx2 * dy2) / factor;
                }
            }
            
            // Apply Neumann BC for pressure correction (∂p'/∂n = 0)
            for i in 0..self.grid.nx {
                workspace.p_prime[(i, 0)] = workspace.p_prime[(i, 1)];
                workspace.p_prime[(i, self.grid.ny-1)] = workspace.p_prime[(i, self.grid.ny-2)];
            }
            for j in 0..self.grid.ny {
                workspace.p_prime[(0, j)] = workspace.p_prime[(1, j)];
                workspace.p_prime[(self.grid.nx-1, j)] = workspace.p_prime[(self.grid.nx-2, j)];
            }
        }
        
        Ok(())
    }

    /// Correct velocity with pressure correction
    fn correct_velocity(
        &self,
        state: &mut SolverState<T>,
        workspace: &PisoWorkspace<T>,
    ) -> Result<()> {
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Velocity correction: u' = -∇p'/a_p
                let grad_p_prime_x = (workspace.p_prime[(i + 1, j)] - workspace.p_prime[(i - 1, j)])
                    / (T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * dx);
                let grad_p_prime_y = (workspace.p_prime[(i, j + 1)] - workspace.p_prime[(i, j - 1)])
                    / (T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? * dy);
                
                state.velocity[(i, j)].x -= grad_p_prime_x / workspace.a_p[(i, j)];
                state.velocity[(i, j)].y -= grad_p_prime_y / workspace.a_p[(i, j)];
                
                // Update pressure
                state.pressure[(i, j)] += workspace.p_prime[(i, j)];
            }
        }
        
        Ok(())
    }

    /// Compute velocity divergence at a cell
    fn compute_velocity_divergence(
        &self,
        state: &SolverState<T>,
        i: usize,
        j: usize,
    ) -> Result<T> {
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let two = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        
        let dudx = (state.velocity[(i + 1, j)].x - state.velocity[(i - 1, j)].x) / (two * dx);
        let dvdy = (state.velocity[(i, j + 1)].y - state.velocity[(i, j - 1)].y) / (two * dy);
        
        Ok(dudx + dvdy)
    }

    /// Apply velocity boundary conditions
    fn apply_velocity_boundary_conditions(&self, state: &mut SolverState<T>) -> Result<()> {
        for (&(i, j), bc) in &self.boundary_conditions {
            match bc {
                BoundaryCondition::Wall { velocity, .. } => {
                    state.velocity[(i, j)] = *velocity;
                }
                BoundaryCondition::VelocityInlet { velocity } => {
                    state.velocity[(i, j)] = *velocity;
                }
                BoundaryCondition::PressureOutlet { .. } => {
                    // Neumann BC for velocity (zero gradient)
                    if i == 0 {
                        state.velocity[(i, j)] = state.velocity[(i + 1, j)];
                    } else if i == self.grid.nx - 1 {
                        state.velocity[(i, j)] = state.velocity[(i - 1, j)];
                    } else if j == 0 {
                        state.velocity[(i, j)] = state.velocity[(i, j + 1)];
                    } else if j == self.grid.ny - 1 {
                        state.velocity[(i, j)] = state.velocity[(i, j - 1)];
                    }
                }
                BoundaryCondition::Symmetry => {
                    // Mirror boundary condition
                    if i == 0 || i == self.grid.nx - 1 {
                        state.velocity[(i, j)].x = T::zero(); // Normal component
                        // Tangential component uses zero gradient
                        let neighbor_i = if i == 0 { i + 1 } else { i - 1 };
                        state.velocity[(i, j)].y = state.velocity[(neighbor_i, j)].y;
                    } else if j == 0 || j == self.grid.ny - 1 {
                        state.velocity[(i, j)].y = T::zero(); // Normal component
                        let neighbor_j = if j == 0 { j + 1 } else { j - 1 };
                        state.velocity[(i, j)].x = state.velocity[(i, neighbor_j)].x;
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Check convergence based on continuity equation
    pub fn check_convergence(&self, state: &SolverState<T>) -> Result<bool> {
        let mut max_continuity_error = T::zero();
        
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                let div_u = self.compute_velocity_divergence(state, i, j)?;
                let error = div_u.abs();
                
                if error > max_continuity_error {
                    max_continuity_error = error;
                }
            }
        }
        
        Ok(max_continuity_error < self.config.pressure_tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_piso_config_builder() {
        let config: PisoConfig<f64> = PisoConfig::builder()
            .time_step(0.001)
            .num_correctors(3)
            .gradient_scheme(GradientScheme::CentralDifference)
            .convection_scheme(ConvectionScheme::QUICK)
            .build()
            .unwrap();
        
        assert_relative_eq!(config.time_step, 0.001, epsilon = 1e-10);
        assert_eq!(config.num_correctors, 3);
        assert!(matches!(config.gradient_scheme, GradientScheme::CentralDifference));
        assert!(matches!(config.convection_scheme, ConvectionScheme::QUICK));
    }

    #[test]
    fn test_solver_state_initialization() {
        let mut state = SolverState::<f64>::new(10, 10);
        state.initialize(Vector2::new(1.0, 0.5), 101325.0);
        
        assert_relative_eq!(state.velocity[(5, 5)].x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(state.velocity[(5, 5)].y, 0.5, epsilon = 1e-10);
        assert_relative_eq!(state.pressure[(5, 5)], 101325.0, epsilon = 1e-10);
    }

    #[test]
    fn test_field_operations() {
        let mut field = ScalarField::<f64>::new(5, 5);
        field.fill(2.5);
        
        assert_relative_eq!(field[(2, 2)], 2.5, epsilon = 1e-10);
        
        field[(2, 2)] = 3.0;
        assert_relative_eq!(field[(2, 2)], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_workspace_creation() {
        let workspace = PisoWorkspace::<f64>::new(10, 10);
        assert_eq!(workspace.p_prime.dimensions(), (10, 10));
        assert_eq!(workspace.a_p.dimensions(), (10, 10));
    }
}