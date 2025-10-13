//! Core momentum equation solver

use super::coefficients::{ConvectionScheme, MomentumCoefficients};
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
    /// Convection discretization scheme
    convection_scheme: ConvectionScheme,
    /// Velocity under-relaxation factor (0 < α ≤ 1, default 0.7)
    /// u_new = α * u_computed + (1-α) * u_old
    velocity_relaxation: T,
}

impl<T: RealField + Copy + FromPrimitive> MomentumSolver<T> {
    /// Create new momentum solver with default deferred correction scheme
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
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(T::one),
        }
    }

    /// Create new momentum solver with specified convection scheme
    pub fn with_convection_scheme(grid: &StructuredGrid2D<T>, scheme: ConvectionScheme) -> Self {
        let config = IterativeSolverConfig::default();
        let linear_solver = BiCGSTAB::new(config);

        Self {
            nx: grid.nx,
            ny: grid.ny,
            dx: grid.dx,
            dy: grid.dy,
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: scheme,
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(T::one),
        }
    }

    /// Set convection scheme
    pub fn set_convection_scheme(&mut self, scheme: ConvectionScheme) {
        self.convection_scheme = scheme;
    }

    /// Set velocity under-relaxation factor
    ///
    /// # Arguments
    /// * `alpha` - Relaxation factor (0 < α ≤ 1)
    ///   - α = 1.0: No relaxation (fastest but may oscillate)
    ///   - α = 0.7: Recommended for most cases
    ///   - α = 0.5: Very stable but slow convergence
    ///
    /// # References
    /// * Patankar (1980) "Numerical Heat Transfer and Fluid Flow" §6.7
    pub fn set_velocity_relaxation(&mut self, alpha: T) {
        self.velocity_relaxation = alpha;
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

        // Diagnostic: Check if coefficients are non-zero (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
            let mut nonzero_count = 0;
            for j in 0..self.ny {
                for i in 0..self.nx {
                    if coeffs.ap.at(i, j).abs() > T::default_epsilon() {
                        nonzero_count += 1;
                    }
                }
            }
            tracing::debug!(
                "Momentum coefficients: {}/{} non-zero entries",
                nonzero_count,
                self.nx * self.ny
            );
        }

        // Assemble linear system
        let (matrix, rhs) = self.assemble_system(&coeffs, component, fields)?;

        // Diagnostic: Check matrix and RHS statistics (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
            let matrix_nnz = matrix.nnz();
            
            tracing::debug!(
                "Linear system: matrix {}x{}, {} nnz",
                matrix.nrows(),
                matrix.ncols(),
                matrix_nnz,
            );
            
            if matrix_nnz == 0 {
                tracing::error!("CRITICAL: Matrix has ZERO non-zero entries - system is empty!");
            }
        }

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
        MomentumCoefficients::compute(
            self.nx,
            self.ny,
            self.dx,
            self.dy,
            dt,
            component,
            fields,
            self.convection_scheme,
        )
    }

    /// Check if a node is on a boundary with Dirichlet BC
    fn is_dirichlet_boundary(&self, i: usize, j: usize) -> bool {
        // Check if node is on any boundary with Dirichlet BC
        if i == 0 && self.boundary_conditions.get("west").map_or(false, |bc| matches!(bc, BoundaryCondition::Dirichlet { .. })) {
            return true;
        }
        if i == self.nx - 1 && self.boundary_conditions.get("east").map_or(false, |bc| matches!(bc, BoundaryCondition::Dirichlet { .. })) {
            return true;
        }
        if j == 0 && self.boundary_conditions.get("south").map_or(false, |bc| matches!(bc, BoundaryCondition::Dirichlet { .. })) {
            return true;
        }
        if j == self.ny - 1 && self.boundary_conditions.get("north").map_or(false, |bc| matches!(bc, BoundaryCondition::Dirichlet { .. })) {
            return true;
        }
        false
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

                // Check if this is a Dirichlet boundary node
                if self.is_dirichlet_boundary(i, j) {
                    // For Dirichlet BC: assemble identity equation φ = bc_value
                    // This is handled in apply_momentum_boundaries, so skip coefficient assembly
                    builder.add_entry(idx, idx, T::one())?;
                    // RHS will be set by boundary condition handler
                    continue;
                }

                // Interior/Neumann nodes: assemble full PDE coefficients
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

        // Apply boundary conditions (sets RHS values for Dirichlet, modifies equations for Neumann)
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
        let alpha = self.velocity_relaxation;
        let one_minus_alpha = T::one() - alpha;

        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;
                let computed_value = solution[idx];

                match component {
                    MomentumComponent::U => {
                        if let Some(u) = fields.u.at_mut(i, j) {
                            let old_value = *u;
                            // Under-relaxation: u_new = α * u_computed + (1-α) * u_old
                            *u = alpha * computed_value + one_minus_alpha * old_value;
                        }
                    }
                    MomentumComponent::V => {
                        if let Some(v) = fields.v.at_mut(i, j) {
                            let old_value = *v;
                            // Under-relaxation: v_new = α * v_computed + (1-α) * v_old
                            *v = alpha * computed_value + one_minus_alpha * old_value;
                        }
                    }
                }
            }
        }
    }
}
