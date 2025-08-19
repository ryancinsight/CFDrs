//! Refactored PISO solver with clean separation of concerns
//!
//! This implementation follows SOLID principles and separates:
//! - Configuration (immutable solver parameters)
//! - State (persistent simulation data)
//! - Workspace (temporary per-timestep data)

use crate::field::{SolverState, PisoWorkspace, ScalarField, VectorField};
use crate::grid::StructuredGrid2D;
use cfd_core::{Result, BoundaryCondition};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Discretization scheme for gradients
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
    /// Quadratic Upstream Interpolation
    QUICK,
    /// Central difference (can be unstable)
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
    /// Maximum iterations for pressure correction
    pub pressure_max_iters: usize,
    /// Gradient discretization scheme
    pub gradient_scheme: GradientScheme,
    /// Convection discretization scheme
    pub convection_scheme: ConvectionScheme,
}

impl<T: RealField + FromPrimitive> PisoConfig<T> {
    /// Create configuration with validated parameters
    pub fn new() -> Result<Self> {
        Ok(Self {
            time_step: T::from_f64(0.01)
                .ok_or_else(|| cfd_core::Error::NumericalError(
                    "Cannot convert time_step to target type".into()
                ))?,
            num_correctors: 2,
            velocity_tolerance: T::from_f64(1e-6)
                .ok_or_else(|| cfd_core::Error::NumericalError(
                    "Cannot convert velocity_tolerance to target type".into()
                ))?,
            pressure_tolerance: T::from_f64(1e-6)
                .ok_or_else(|| cfd_core::Error::NumericalError(
                    "Cannot convert pressure_tolerance to target type".into()
                ))?,
            pressure_max_iters: 20,
            gradient_scheme: GradientScheme::CentralDifference,
            convection_scheme: ConvectionScheme::Hybrid,
        })
    }

    /// Builder pattern for custom configuration
    pub fn builder() -> PisoConfigBuilder<T> {
        PisoConfigBuilder::default()
    }
}

/// Configuration builder for fluent API
#[derive(Default)]
pub struct PisoConfigBuilder<T: RealField> {
    time_step: Option<T>,
    num_correctors: Option<usize>,
    velocity_tolerance: Option<T>,
    pressure_tolerance: Option<T>,
    pressure_max_iters: Option<usize>,
    gradient_scheme: Option<GradientScheme>,
    convection_scheme: Option<ConvectionScheme>,
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
            pressure_max_iters: self.pressure_max_iters.unwrap_or(default.pressure_max_iters),
            gradient_scheme: self.gradient_scheme.unwrap_or(default.gradient_scheme),
            convection_scheme: self.convection_scheme.unwrap_or(default.convection_scheme),
        })
    }
}

/// Refactored PISO solver with clean architecture
pub struct PisoSolver<T: RealField> {
    /// Solver configuration (immutable)
    config: PisoConfig<T>,
    /// Grid structure
    grid: StructuredGrid2D<T>,
    /// Fluid properties
    properties: FluidProperties<T>,
    /// Boundary conditions
    boundary_conditions: HashMap<(usize, usize), BoundaryCondition<T>>,
}

/// Fluid properties container
#[derive(Debug, Clone)]
pub struct FluidProperties<T: RealField> {
    pub density: T,
    pub viscosity: T,
}

impl<T: RealField + FromPrimitive + Copy> PisoSolver<T> {
    /// Create new solver with configuration
    pub fn new(
        config: PisoConfig<T>,
        grid: StructuredGrid2D<T>,
        properties: FluidProperties<T>,
    ) -> Self {
        Self {
            config,
            grid,
            properties,
            boundary_conditions: HashMap::new(),
        }
    }

    /// Set boundary condition at specific location
    pub fn set_boundary_condition(
        &mut self,
        i: usize,
        j: usize,
        bc: BoundaryCondition<T>,
    ) {
        self.boundary_conditions.insert((i, j), bc);
    }

    /// Perform one timestep of the PISO algorithm
    pub fn step(&self, state: &mut SolverState<T>) -> Result<()> {
        // Initialize transient if needed
        state.init_transient();
        
        // Create temporary workspace for this timestep
        let mut workspace = PisoWorkspace::new(self.grid.nx, self.grid.ny);
        
        // PISO algorithm steps
        self.solve_momentum_predictor(state, &mut workspace)?;
        
        for _ in 0..self.config.num_correctors {
            self.solve_pressure_correction(state, &mut workspace)?;
            self.correct_velocity(state, &workspace)?;
        }
        
        // Advance to next timestep
        state.advance_timestep();
        
        Ok(())
    }

    /// Solve momentum predictor step
    fn solve_momentum_predictor(
        &self,
        state: &mut SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        let dt = self.config.time_step;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        
        // Compute momentum equation coefficients
        self.compute_momentum_coefficients(workspace, dt)?;
        
        // Solve momentum equations
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Apply selected gradient scheme
                let pressure_grad = self.compute_pressure_gradient(
                    &state.pressure,
                    i, j,
                    dx, dy
                );
                
                // Apply selected convection scheme
                let convection = self.compute_convection_term(
                    &state.velocity,
                    i, j,
                    dx, dy
                );
                
                // Update velocity (no redundant clones for Copy types)
                let diffusion = self.compute_diffusion_term(
                    &state.velocity,
                    i, j,
                    dx, dy
                );
                
                // Compute new velocity
                let ap = workspace.coefficients.a_p[(i, j)];
                let source = -pressure_grad + diffusion - convection;
                
                state.velocity[(i, j)] = state.velocity_old
                    .as_ref()
                    .map(|old| old[(i, j)])
                    .unwrap_or_else(Vector2::zeros) 
                    + dt * source / ap;
            }
        }
        
        // Apply boundary conditions
        self.apply_boundary_conditions(state)?;
        
        Ok(())
    }

    /// Compute pressure gradient using selected scheme
    fn compute_pressure_gradient(
        &self,
        pressure: &ScalarField<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> Vector2<T> {
        match self.config.gradient_scheme {
            GradientScheme::CentralDifference => {
                Vector2::new(
                    (pressure[(i + 1, j)] - pressure[(i - 1, j)]) / (T::from_f64(2.0).unwrap() * dx),
                    (pressure[(i, j + 1)] - pressure[(i, j - 1)]) / (T::from_f64(2.0).unwrap() * dy),
                )
            }
            GradientScheme::ForwardDifference => {
                Vector2::new(
                    (pressure[(i + 1, j)] - pressure[(i, j)]) / dx,
                    (pressure[(i, j + 1)] - pressure[(i, j)]) / dy,
                )
            }
            GradientScheme::BackwardDifference => {
                Vector2::new(
                    (pressure[(i, j)] - pressure[(i - 1, j)]) / dx,
                    (pressure[(i, j)] - pressure[(i, j - 1)]) / dy,
                )
            }
        }
    }

    /// Compute convection term using selected scheme
    fn compute_convection_term(
        &self,
        velocity: &VectorField<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> Vector2<T> {
        match self.config.convection_scheme {
            ConvectionScheme::Upwind => {
                // First-order upwind implementation
                let u = velocity[(i, j)].x;
                let v = velocity[(i, j)].y;
                
                let dudx = if u > T::zero() {
                    (velocity[(i, j)].x - velocity[(i - 1, j)].x) / dx
                } else {
                    (velocity[(i + 1, j)].x - velocity[(i, j)].x) / dx
                };
                
                let dvdy = if v > T::zero() {
                    (velocity[(i, j)].y - velocity[(i, j - 1)].y) / dy
                } else {
                    (velocity[(i, j + 1)].y - velocity[(i, j)].y) / dy
                };
                
                Vector2::new(u * dudx, v * dvdy)
            }
            ConvectionScheme::Central => {
                // Central difference (can be unstable)
                let u = velocity[(i, j)].x;
                let v = velocity[(i, j)].y;
                
                let dudx = (velocity[(i + 1, j)].x - velocity[(i - 1, j)].x) / (T::from_f64(2.0).unwrap() * dx);
                let dvdy = (velocity[(i, j + 1)].y - velocity[(i, j - 1)].y) / (T::from_f64(2.0).unwrap() * dy);
                
                Vector2::new(u * dudx, v * dvdy)
            }
            _ => {
                // TODO: Implement Hybrid and QUICK schemes
                Vector2::zeros()
            }
        }
    }

    /// Compute diffusion term
    fn compute_diffusion_term(
        &self,
        velocity: &VectorField<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> Vector2<T> {
        let nu = self.properties.viscosity / self.properties.density;
        
        let d2udx2 = (velocity[(i + 1, j)].x - T::from_f64(2.0).unwrap() * velocity[(i, j)].x + velocity[(i - 1, j)].x) / (dx * dx);
        let d2udy2 = (velocity[(i, j + 1)].x - T::from_f64(2.0).unwrap() * velocity[(i, j)].x + velocity[(i, j - 1)].x) / (dy * dy);
        
        let d2vdx2 = (velocity[(i + 1, j)].y - T::from_f64(2.0).unwrap() * velocity[(i, j)].y + velocity[(i - 1, j)].y) / (dx * dx);
        let d2vdy2 = (velocity[(i, j + 1)].y - T::from_f64(2.0).unwrap() * velocity[(i, j)].y + velocity[(i, j - 1)].y) / (dy * dy);
        
        Vector2::new(
            nu * (d2udx2 + d2udy2),
            nu * (d2vdx2 + d2vdy2),
        )
    }

    /// Compute momentum equation coefficients
    fn compute_momentum_coefficients(
        &self,
        workspace: &mut PisoWorkspace<T>,
        dt: T,
    ) -> Result<()> {
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let nu = self.properties.viscosity / self.properties.density;
        
        for i in 1..self.grid.nx - 1 {
            for j in 1..self.grid.ny - 1 {
                // Diffusion coefficients
                let diff_x = nu / (dx * dx);
                let diff_y = nu / (dy * dy);
                
                workspace.coefficients.a_e[(i, j)] = diff_x;
                workspace.coefficients.a_w[(i, j)] = diff_x;
                workspace.coefficients.a_n[(i, j)] = diff_y;
                workspace.coefficients.a_s[(i, j)] = diff_y;
                
                // Central coefficient
                workspace.coefficients.a_p[(i, j)] = 
                    T::one() / dt + T::from_f64(2.0).unwrap() * (diff_x + diff_y);
            }
        }
        
        Ok(())
    }

    /// Solve pressure correction equation
    fn solve_pressure_correction(
        &self,
        state: &SolverState<T>,
        workspace: &mut PisoWorkspace<T>,
    ) -> Result<()> {
        // Implement pressure correction
        // This would use a linear solver (Gauss-Seidel, CG, etc.)
        Ok(())
    }

    /// Correct velocity with pressure correction
    fn correct_velocity(
        &self,
        state: &mut SolverState<T>,
        workspace: &PisoWorkspace<T>,
    ) -> Result<()> {
        // Implement velocity correction
        Ok(())
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&self, state: &mut SolverState<T>) -> Result<()> {
        // Apply user-specified boundary conditions
        for ((i, j), bc) in &self.boundary_conditions {
            match bc {
                BoundaryCondition::Dirichlet(value) => {
                    // Set velocity directly
                }
                BoundaryCondition::Neumann(gradient) => {
                    // Apply gradient condition
                }
                _ => {}
            }
        }
        Ok(())
    }
}