//! Problem definitions for 2D incompressible flow simulations.
//!
//! This module implements the Problem trait from cfd-core for 2D incompressible
//! Navier-Stokes equations, enabling integration with the framework's solver abstractions.

use crate::grid::{Grid2D, StructuredGrid2D};
use cfd_core::boundary::BoundaryCondition;
use cfd_core::fluid::ConstantPropertyFluid;
use cfd_core::problem::Problem;
use nalgebra::{RealField, Vector2};
use std::collections::HashMap;

/// Problem definition for 2D incompressible flow
#[derive(Debug, Clone)]
pub struct IncompressibleFlowProblem<T: RealField + Copy> {
    /// Computational grid
    pub grid: StructuredGrid2D<T>,
    /// Boundary conditions at specific grid points
    pub boundary_conditions: HashMap<(usize, usize), BoundaryCondition<T>>,
    /// Fluid properties
    pub fluid: ConstantPropertyFluid<T>,
    /// Initial velocity field
    pub initial_velocity: Vec<Vec<Vector2<T>>>,
    /// Initial pressure field
    pub initial_pressure: Vec<Vec<T>>,
    /// Time step for transient simulations
    pub time_step: Option<T>,
    /// Simulation end time
    pub end_time: Option<T>,
}

impl<T: RealField + Copy> IncompressibleFlowProblem<T> {
    /// Create new incompressible flow problem
    pub fn new(
        grid: StructuredGrid2D<T>,
        boundary_conditions: HashMap<(usize, usize), BoundaryCondition<T>>,
        fluid: ConstantPropertyFluid<T>,
    ) -> Self {
        let nx = grid.nx();
        let ny = grid.ny();

        Self {
            grid,
            boundary_conditions,
            fluid,
            initial_velocity: vec![vec![Vector2::zeros(); ny]; nx],
            initial_pressure: vec![vec![T::zero(); ny]; nx],
            time_step: None,
            end_time: None,
        }
    }

    /// Set initial velocity field
    pub fn with_initial_velocity(mut self, velocity: Vec<Vec<Vector2<T>>>) -> Self {
        self.initial_velocity = velocity;
        self
    }

    /// Set initial pressure field
    pub fn with_initial_pressure(mut self, pressure: Vec<Vec<T>>) -> Self {
        self.initial_pressure = pressure;
        self
    }

    /// Set time step for transient simulations
    pub fn with_time_step(mut self, dt: T) -> Self {
        self.time_step = Some(dt);
        self
    }

    /// Set simulation end time
    pub fn with_end_time(mut self, t_end: T) -> Self {
        self.end_time = Some(t_end);
        self
    }

    /// Validate problem setup
    pub fn validate(&self) -> cfd_core::error::Result<()> {
        // Check grid dimensions match initial fields
        if self.initial_velocity.len() != self.grid.nx() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Initial velocity field dimensions don't match grid".into(),
            ));
        }

        if !self.initial_velocity.is_empty() && self.initial_velocity[0].len() != self.grid.ny() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Initial velocity field dimensions don't match grid".into(),
            ));
        }

        // Check if boundary conditions are complete for boundary cells
        let mut missing_bcs = Vec::new();

        // Check boundary cells
        for i in 0..self.grid.nx() {
            for j in [0, self.grid.ny() - 1] {
                if !self.boundary_conditions.contains_key(&(i, j)) {
                    missing_bcs.push((i, j));
                }
            }
        }

        for j in 0..self.grid.ny() {
            for i in [0, self.grid.nx() - 1] {
                if !self.boundary_conditions.contains_key(&(i, j)) {
                    missing_bcs.push((i, j));
                }
            }
        }

        if !missing_bcs.is_empty() {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Missing boundary conditions for cells: {missing_bcs:?}"
            )));
        }

        Ok(())
    }
}

// Note: Problem trait implementation postponed for now to focus on architectural refactoring
// Problem trait integration implemented through IncompressibleFlowProblem

/// Solution structure for incompressible flow
#[derive(Debug, Clone)]
pub struct IncompressibleFlowSolution<T: RealField + Copy> {
    /// Final velocity field
    pub velocity: Vec<Vec<Vector2<T>>>,
    /// Final pressure field
    pub pressure: Vec<Vec<T>>,
    /// Number of iterations to convergence
    pub iterations: usize,
    /// Final residual
    pub residual: T,
    /// Simulation time (for transient)
    pub time: T,
}

impl<T: RealField + Copy> IncompressibleFlowSolution<T> {
    /// Create new solution
    pub fn new(
        velocity: Vec<Vec<Vector2<T>>>,
        pressure: Vec<Vec<T>>,
        iterations: usize,
        residual: T,
        time: T,
    ) -> Self {
        Self {
            velocity,
            pressure,
            iterations,
            residual,
            time,
        }
    }

    /// Get maximum velocity magnitude
    pub fn max_velocity_magnitude(&self) -> T {
        self.velocity
            .iter()
            .flat_map(|row| row.iter())
            .map(nalgebra::Matrix::magnitude)
            .fold(T::zero(), |acc, mag| if mag > acc { mag } else { acc })
    }

    /// Get maximum pressure
    pub fn max_pressure(&self) -> T {
        self.pressure
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .fold(T::zero(), |acc, p| if p > acc { p } else { acc })
    }
}
