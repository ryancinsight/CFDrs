//! Core momentum equation solver

use super::coefficients::MomentumCoefficients;
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use cfd_core::boundary::BoundaryCondition;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::IterativeSolverConfig;
use cfd_math::linear_solver::{BiCGSTAB, IterativeLinearSolver};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Component of momentum equation (U or V)
#[derive(Debug, Clone, Copy)]
pub enum MomentumComponent {
    /// U-component (x-direction velocity)
    U,
    /// V-component (y-direction velocity)
    V,
}

/// Momentum equation solver for 2D incompressible flow
pub struct MomentumSolver<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Boundary conditions
    boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Linear solver
    linear_solver: BiCGSTAB<T>,
}

impl<T: RealField + Copy + FromPrimitive> MomentumSolver<T> {
    /// Create new momentum solver
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        let config = IterativeSolverConfig::default();
        let linear_solver = BiCGSTAB::new(config);

        Self {
            nx: grid.nx,
            ny: grid.ny,
            dx: grid.dx,
            dy: grid.dy,
            boundary_conditions: HashMap::new(),
            linear_solver,
        }
    }

    /// Set boundary condition
    pub fn set_boundary(&mut self, name: String, bc: BoundaryCondition<T>) {
        self.boundary_conditions.insert(name, bc);
    }

    /// Solve momentum equation for specified component
    pub fn solve(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        // Compute coefficients
        let coeffs = self.compute_coefficients(component, fields, dt)?;

        // Assemble linear system
        let (matrix, rhs) = self.assemble_system(&coeffs, component, fields)?;

        // Solve linear system
        let mut solution = DVector::zeros(matrix.nrows());
        self.linear_solver.solve(
            &matrix,
            &rhs,
            &mut solution,
            None::<&IdentityPreconditioner>,
        )?;

        // Update velocity field
        self.update_velocity(component, fields, &solution);

        Ok(())
    }

    fn compute_coefficients(
        &self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<MomentumCoefficients<T>> {
        MomentumCoefficients::compute(self.nx, self.ny, self.dx, self.dy, dt, component, fields)
    }

    fn assemble_system(
        &self,
        coeffs: &MomentumCoefficients<T>,
        component: MomentumComponent,
        _fields: &SimulationFields<T>,
    ) -> cfd_core::error::Result<(SparseMatrix<T>, DVector<T>)> {
        let n = self.nx * self.ny;
        let mut builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::zeros(n);

        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;

                // Central coefficient
                builder.add_entry(idx, idx, coeffs.ap.at(i, j))?;

                // Neighbor coefficients
                if i > 0 {
                    builder.add_entry(idx, idx - 1, -coeffs.aw.at(i, j))?;
                }
                if i < self.nx - 1 {
                    builder.add_entry(idx, idx + 1, -coeffs.ae.at(i, j))?;
                }
                if j > 0 {
                    builder.add_entry(idx, idx - self.nx, -coeffs.as_.at(i, j))?;
                }
                if j < self.ny - 1 {
                    builder.add_entry(idx, idx + self.nx, -coeffs.an.at(i, j))?;
                }

                // Source term including pressure gradient
                rhs[idx] = coeffs.source.at(i, j);
            }
        }

        // Apply boundary conditions
        super::boundary::apply_momentum_boundaries(
            &mut builder,
            &mut rhs,
            component,
            &self.boundary_conditions,
            self.nx,
            self.ny,
        )?;

        Ok((builder.build()?, rhs))
    }

    fn update_velocity(
        &self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        solution: &DVector<T>,
    ) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;
                match component {
                    MomentumComponent::U => {
                        if let Some(u) = fields.u.at_mut(i, j) {
                            *u = solution[idx];
                        }
                    }
                    MomentumComponent::V => {
                        if let Some(v) = fields.v.at_mut(i, j) {
                            *v = solution[idx];
                        }
                    }
                }
            }
        }
    }
}
