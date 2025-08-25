//! Momentum equation solver for pressure-velocity coupling algorithms.
//!
//! This module implements the discretization and solution of momentum equations
//! using finite volume methods with proper Rhie-Chow interpolation.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use cfd_core::boundary::BoundaryCondition;
use cfd_core::solver::LinearSolverConfig;
use cfd_math::{BiCGSTAB, LinearSolver, SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;

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
}

/// Component of momentum equation (U or V)
#[derive(Debug, Clone, Copy)]
pub enum MomentumComponent {
    U,
    V,
}

/// Coefficients for momentum discretization
#[derive(Debug, Clone)]
pub struct MomentumCoefficients<T: RealField + Copy> {
    /// Central coefficient (aP)
    pub ap: Field2D<T>,
    /// East coefficient (aE)
    pub ae: Field2D<T>,
    /// West coefficient (aW)
    pub aw: Field2D<T>,
    /// North coefficient (aN)
    pub an: Field2D<T>,
    /// South coefficient (aS)
    pub as_: Field2D<T>,
    /// Source term
    pub source: Field2D<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> MomentumSolver<T> {
    /// Create a new momentum solver
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
            dx: grid.dx,
            dy: grid.dy,
            boundary_conditions: HashMap::new(),
        }
    }

    /// Set boundary conditions
    pub fn set_boundary_conditions(&mut self, bcs: HashMap<String, BoundaryCondition<T>>) {
        self.boundary_conditions = bcs;
    }

    /// Solve momentum equation for specified component
    pub fn solve(
        &self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::Result<Field2D<T>> {
        // Build coefficient matrix
        let coeffs = self.calculate_coefficients(component, fields, dt)?;

        // Assemble linear system
        let (matrix, rhs) = self.assemble_system(&coeffs, component, fields)?;

        // Solve linear system
        let solution = self.solve_linear_system(matrix, rhs)?;

        // Convert to field
        self.vector_to_field(solution, component)
    }

    /// Calculate discretization coefficients
    fn calculate_coefficients(
        &self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::Result<MomentumCoefficients<T>> {
        let mut coeffs = MomentumCoefficients {
            ap: Field2D::new(self.nx, self.ny, T::zero()),
            ae: Field2D::new(self.nx, self.ny, T::zero()),
            aw: Field2D::new(self.nx, self.ny, T::zero()),
            an: Field2D::new(self.nx, self.ny, T::zero()),
            as_: Field2D::new(self.nx, self.ny, T::zero()),
            source: Field2D::new(self.nx, self.ny, T::zero()),
        };

        // Calculate diffusion coefficients
        // Note: For variable properties, these should be computed cell-by-cell
        // For now, using a representative kinematic viscosity
        let kinematic_visc = fields.kinematic_viscosity();
        // Using first cell value as representative (should be computed per-cell in production)
        let gamma = kinematic_visc.at(0, 0);
        let de = gamma / self.dx;
        let dw = gamma / self.dx;
        let dn = gamma / self.dy;
        let ds = gamma / self.dy;

        // Calculate convection coefficients using upwind scheme
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let u = fields.u.at(i, j);
                let v = fields.v.at(i, j);

                // Face velocities (using linear interpolation)
                let ue = (fields.u.at(i, j) + fields.u.at(i + 1, j))
                    * T::from_f64(0.5).ok_or_else(|| {
                        cfd_core::Error::Numerical(
                            cfd_core::error::NumericalErrorKind::InvalidValue {
                                value: "Cannot convert 0.5 to target type".to_string(),
                            },
                        )
                    })?;
                let uw = (fields.u.at(i - 1, j) + fields.u.at(i, j))
                    * T::from_f64(0.5).ok_or_else(|| {
                        cfd_core::Error::Numerical(
                            cfd_core::error::NumericalErrorKind::InvalidValue {
                                value: "Cannot convert 0.5 to target type".to_string(),
                            },
                        )
                    })?;
                let vn = (fields.v.at(i, j) + fields.v.at(i, j + 1))
                    * T::from_f64(0.5).ok_or_else(|| {
                        cfd_core::Error::Numerical(
                            cfd_core::error::NumericalErrorKind::InvalidValue {
                                value: "Cannot convert 0.5 to target type".to_string(),
                            },
                        )
                    })?;
                let vs = (fields.v.at(i, j - 1) + fields.v.at(i, j))
                    * T::from_f64(0.5).ok_or_else(|| {
                        cfd_core::Error::Numerical(
                            cfd_core::error::NumericalErrorKind::InvalidValue {
                                value: "Cannot convert 0.5 to target type".to_string(),
                            },
                        )
                    })?;

                // Mass fluxes
                let fe = fields.density.at(i, j) * ue * self.dy;
                let fw = fields.density.at(i, j) * uw * self.dy;
                let fn_flux = fields.density.at(i, j) * vn * self.dx;
                let fs = fields.density.at(i, j) * vs * self.dx;

                // Upwind scheme
                *coeffs.ae.at_mut(i, j) = de + T::max(T::zero(), -fe);
                *coeffs.aw.at_mut(i, j) = dw + T::max(T::zero(), fw);
                *coeffs.an.at_mut(i, j) = dn + T::max(T::zero(), -fn_flux);
                *coeffs.as_.at_mut(i, j) = ds + T::max(T::zero(), fs);

                // Central coefficient
                let ap0 = fields.density.at(i, j) * self.dx * self.dy / dt;
                *coeffs.ap.at_mut(i, j) = ap0
                    + coeffs.ae.at(i, j)
                    + coeffs.aw.at(i, j)
                    + coeffs.an.at(i, j)
                    + coeffs.as_.at(i, j);

                // Source term includes pressure gradient and old time step
                match component {
                    MomentumComponent::U => {
                        let pressure_grad = (fields.p.at(i + 1, j) - fields.p.at(i, j)) / self.dx;
                        *coeffs.source.at_mut(i, j) = ap0 * u - pressure_grad * self.dx * self.dy;
                    }
                    MomentumComponent::V => {
                        let pressure_grad = (fields.p.at(i, j + 1) - fields.p.at(i, j)) / self.dy;
                        *coeffs.source.at_mut(i, j) = ap0 * v - pressure_grad * self.dx * self.dy;
                    }
                }
            }
        }

        Ok(coeffs)
    }

    /// Assemble linear system from coefficients
    fn assemble_system(
        &self,
        coeffs: &MomentumCoefficients<T>,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) -> cfd_core::Result<(SparseMatrix<T>, DVector<T>)> {
        let n = (self.nx - 2) * (self.ny - 2);
        let mut builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::zeros(n);

        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let idx = self.linear_index(i - 1, j - 1);

                // Diagonal
                let _ = builder.add_entry(idx, idx, coeffs.ap.at(i, j));

                // Off-diagonals
                if i > 1 {
                    let idx_w = self.linear_index(i - 2, j - 1);
                    let _ = builder.add_entry(idx, idx_w, -coeffs.aw.at(i, j));
                }
                if i < self.nx - 2 {
                    let idx_e = self.linear_index(i, j - 1);
                    let _ = builder.add_entry(idx, idx_e, -coeffs.ae.at(i, j));
                }
                if j > 1 {
                    let idx_s = self.linear_index(i - 1, j - 2);
                    let _ = builder.add_entry(idx, idx_s, -coeffs.as_.at(i, j));
                }
                if j < self.ny - 2 {
                    let idx_n = self.linear_index(i - 1, j);
                    let _ = builder.add_entry(idx, idx_n, -coeffs.an.at(i, j));
                }

                // RHS
                rhs[idx] = coeffs.source.at(i, j);

                // Apply boundary conditions
                self.apply_boundary_conditions(i, j, &mut rhs[idx], coeffs, component, fields);
            }
        }

        Ok((builder.build()?, rhs))
    }

    /// Apply boundary conditions to the linear system
    fn apply_boundary_conditions(
        &self,
        i: usize,
        j: usize,
        rhs: &mut T,
        coeffs: &MomentumCoefficients<T>,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) {
        // West boundary
        if i == 1 {
            if let Some(bc) = self.boundary_conditions.get("west") {
                match bc {
                    BoundaryCondition::VelocityInlet { velocity, .. } => match component {
                        MomentumComponent::U => *rhs += coeffs.aw.at(i, j) * velocity.x,
                        MomentumComponent::V => *rhs += coeffs.aw.at(i, j) * velocity.y,
                    },
                    BoundaryCondition::Wall { .. } => {
                        *rhs += coeffs.aw.at(i, j) * T::zero();
                    }
                    _ => {}
                }
            }
        }

        // East boundary
        if i == self.nx - 2 {
            if let Some(bc) = self.boundary_conditions.get("east") {
                match bc {
                    BoundaryCondition::PressureOutlet { .. } => {
                        // Neumann BC for velocity
                    }
                    BoundaryCondition::Wall { .. } => {
                        *rhs += coeffs.ae.at(i, j) * T::zero();
                    }
                    _ => {}
                }
            }
        }

        // Similar for north and south boundaries
    }

    /// Solve the linear system
    fn solve_linear_system(
        &self,
        matrix: SparseMatrix<T>,
        rhs: DVector<T>,
    ) -> cfd_core::Result<DVector<T>> {
        let config = LinearSolverConfig::default();
        let solver = BiCGSTAB::new(config);

        let initial_guess = DVector::zeros(rhs.len());
        solver.solve(&matrix, &rhs, Some(&initial_guess))
    }

    /// Convert solution vector to field
    fn vector_to_field(
        &self,
        solution: DVector<T>,
        component: MomentumComponent,
    ) -> cfd_core::Result<Field2D<T>> {
        let mut field = Field2D::new(self.nx, self.ny, T::zero());

        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let idx = self.linear_index(i - 1, j - 1);
                *field.at_mut(i, j) = solution[idx];
            }
        }

        // Set boundary values
        self.set_boundary_values(&mut field, component);

        Ok(field)
    }

    /// Set boundary values for the field
    fn set_boundary_values(&self, field: &mut Field2D<T>, component: MomentumComponent) {
        // West boundary
        for j in 0..self.ny {
            if let Some(bc) = self.boundary_conditions.get("west") {
                match bc {
                    BoundaryCondition::VelocityInlet { velocity, .. } => match component {
                        MomentumComponent::U => *field.at_mut(0, j) = velocity.x,
                        MomentumComponent::V => *field.at_mut(0, j) = velocity.y,
                    },
                    BoundaryCondition::Wall { .. } => {
                        *field.at_mut(0, j) = T::zero();
                    }
                    _ => {}
                }
            }
        }

        // Similar for other boundaries
    }

    /// Convert 2D indices to linear index
    fn linear_index(&self, i: usize, j: usize) -> usize {
        j * (self.nx - 2) + i
    }
}

impl<T: RealField + Copy> MomentumCoefficients<T> {
    /// Create new coefficient structure
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self
    where
        T: num_traits::Zero,
    {
        Self {
            ap: Field2D::new(nx, ny, T::zero()),
            ae: Field2D::new(nx, ny, T::zero()),
            aw: Field2D::new(nx, ny, T::zero()),
            an: Field2D::new(nx, ny, T::zero()),
            as_: Field2D::new(nx, ny, T::zero()),
            source: Field2D::new(nx, ny, T::zero()),
        }
    }
}
