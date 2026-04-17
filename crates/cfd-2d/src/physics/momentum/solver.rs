//! Core momentum equation solver
//!
//! # Theorem
//! The momentum discretization must conserve linear momentum globally and locally.
//!
//! **Proof sketch**:
//! By integrating the Navier-Stokes momentum equation over a control volume $\Omega$,
//! Gauss's divergence theorem converts the convective and diffusive volume integrals
//! into surface fluxes. The finite volume method ensures that the flux leaving one
//! cell exactly equals the flux entering the adjacent cell. Thus, in the absence of
//! external forces and boundary fluxes, the total momentum $\int_\Omega \rho \mathbf{u} dV$
//! is exactly conserved to machine precision.

use super::coefficients::{ConvectionScheme, MomentumCoefficients};
use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::turbulence::TurbulenceModel;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::DirectSparseSolver;
use cfd_math::linear_solver::IterativeSolverConfig;
use cfd_math::linear_solver::{IterativeLinearSolver, GMRES};

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
    /// Grid reference for boundary condition calculations
    grid: StructuredGrid2D<T>,
    /// Boundary conditions
    boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Linear solver
    linear_solver: GMRES<T>,

    /// Convection discretization scheme
    convection_scheme: ConvectionScheme,
    /// Velocity under-relaxation factor (0 < α ≤ 1, default 0.7)
    /// u_new = α * u_computed + (1-α) * u_old
    velocity_relaxation: T,
    /// Optional turbulence model for computing turbulent viscosity
    turbulence_model: Option<Box<dyn TurbulenceModel<T>>>,

    /// U-momentum coefficients
    coeffs_u: MomentumCoefficients<T>,
    /// V-momentum coefficients
    coeffs_v: MomentumCoefficients<T>,
    
    matrix_u: Option<SparseMatrix<T>>,
    matrix_v: Option<SparseMatrix<T>>,
    matrix_builder_u: Option<SparseMatrixBuilder<T>>,
    matrix_builder_v: Option<SparseMatrixBuilder<T>>,
    rhs_u: Option<DVector<T>>,
    rhs_v: Option<DVector<T>>,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::ToPrimitive> MomentumSolver<T> {
    /// Create new momentum solver with default deferred correction scheme
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        let linear_solver = GMRES::new(Self::linear_solver_config(), 30);

        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: T::from_f64(0.7).expect("analytical constant conversion"),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::U, &SimulationFields::new(grid.nx, grid.ny), ConvectionScheme::default()).unwrap(),
            coeffs_v: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::V, &SimulationFields::new(grid.nx, grid.ny), ConvectionScheme::default()).unwrap(),
            matrix_u: None,
            matrix_v: None,
            matrix_builder_u: None,
            matrix_builder_v: None,
            rhs_u: None,
            rhs_v: None,
        }
    }

    /// Create new momentum solver with specified convection scheme
    pub fn with_convection_scheme(grid: &StructuredGrid2D<T>, scheme: ConvectionScheme) -> Self {
        let linear_solver = GMRES::new(Self::linear_solver_config(), 30);

        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: scheme,
            velocity_relaxation: T::from_f64(0.7).expect("analytical constant conversion"),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::U, &SimulationFields::new(grid.nx, grid.ny), scheme).unwrap(),
            coeffs_v: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::V, &SimulationFields::new(grid.nx, grid.ny), scheme).unwrap(),
            matrix_u: None,
            matrix_v: None,
            matrix_builder_u: None,
            matrix_builder_v: None,
            rhs_u: None,
            rhs_v: None,
        }
    }

    /// Create new momentum solver with parallel SpMV enabled for multi-core performance
    #[must_use]
    pub fn with_parallel_spmv(grid: &StructuredGrid2D<T>) -> Self {
        let mut config = Self::linear_solver_config();
        config.use_parallel_spmv = true;
        let linear_solver = GMRES::new(config, 30);

        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: T::from_f64(0.7).expect("analytical constant conversion"),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::U, &SimulationFields::new(grid.nx, grid.ny), ConvectionScheme::default()).unwrap(),
            coeffs_v: MomentumCoefficients::compute(grid.nx, grid.ny, T::one(), T::one(), T::one(), MomentumComponent::V, &SimulationFields::new(grid.nx, grid.ny), ConvectionScheme::default()).unwrap(),
            matrix_u: None,
            matrix_v: None,
            matrix_builder_u: None,
            matrix_builder_v: None,
            rhs_u: None,
            rhs_v: None,
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

    /// Get boundary conditions
    #[must_use]
    pub fn boundary_conditions(&self) -> &HashMap<String, BoundaryCondition<T>> {
        &self.boundary_conditions
    }

    /// Validate boundary condition consistency
    pub fn validate_boundary_conditions(&self) -> cfd_core::error::Result<()> {
        super::boundary::validate_boundary_consistency(&self.boundary_conditions, &self.grid)
            .map_err(|e| cfd_core::error::Error::InvalidConfiguration(e.to_string()))
    }

    /// Set turbulence model for computing turbulent viscosity
    pub fn set_turbulence_model(&mut self, model: Box<dyn TurbulenceModel<T>>) {
        self.turbulence_model = Some(model);
    }

    /// Remove turbulence model (use laminar flow only)
    pub fn remove_turbulence_model(&mut self) {
        self.turbulence_model = None;
    }

    /// Check if turbulence model is active
    #[must_use]
    pub fn has_turbulence_model(&self) -> bool {
        self.turbulence_model.is_some()
    }

    /// Solve momentum equation for specified component.
    ///
    /// Integrates with any configured turbulence model by updating
    /// estimated turbulence quantities based on flow physics.
    pub fn solve(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        self.solve_with_coefficients(component, fields, dt)?;
        Ok(())
    }

    /// Solve momentum equation and return coefficients for Rhie-Chow interpolation
    pub fn solve_with_coefficients(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        // Compute coefficients
        self.compute_coefficients_into(component, fields, dt)?;
        
        // Diagnostic: Check if coefficients are non-zero (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
            let coeffs = match component {
                MomentumComponent::U => &self.coeffs_u,
                MomentumComponent::V => &self.coeffs_v,
            };
            let mut nonzero_count = 0;
            for j in 0..self.grid.ny {
                for i in 0..self.grid.nx {
                    if coeffs.ap.at(i, j).abs() > T::default_epsilon() {
                        nonzero_count += 1;
                    }
                }
            }
            tracing::debug!(
                "Momentum coefficients: {}/{} non-zero entries",
                nonzero_count,
                self.grid.nx * self.grid.ny
            );
        }

        // Assemble linear system topology once and update values
        self.assemble_system_into(component, fields)?;
        
        let (matrix, rhs) = match component {
            MomentumComponent::U => (self.matrix_u.as_ref().unwrap(), self.rhs_u.as_ref().unwrap()),
            MomentumComponent::V => (self.matrix_v.as_ref().unwrap(), self.rhs_v.as_ref().unwrap()),
        };

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
        let _rhs_norm = rhs.norm();
        let solve_result = self.linear_solver.solve(
            matrix,
            rhs,
            &mut solution,
            None::<&IdentityPreconditioner>,
        );

        match solve_result {
            Ok(_) => {}
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { max },
            )) => {
                tracing::warn!(
                    component = ?component,
                    max_iterations = max,
                    "Momentum GMRES stalled; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Momentum direct sparse solve failed for {component:?}: {error}"
                    ))
                })?;
            }
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Breakdown,
            )) => {
                tracing::warn!(
                    component = ?component,
                    "Momentum GMRES breakdown; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Momentum direct sparse solve failed for {component:?}: {error}"
                    ))
                })?;
            }
            Err(error) => return Err(error),
        }

        // Update velocity field
        self.update_velocity(component, fields, &solution);


        Ok(())
    }

    fn linear_solver_config() -> IterativeSolverConfig<T> {
        IterativeSolverConfig {
            max_iterations: crate::constants::solver::DEFAULT_MAX_ITERATIONS,
            tolerance: T::from_f64(crate::constants::solver::DEFAULT_TOLERANCE)
                .expect("Failed to convert momentum solver tolerance"),
            use_preconditioner: false,
            use_parallel_spmv: false,
        }
    }

    /// Get the last computed A_p and A_p_consistent coefficients for both velocity components
    pub fn get_ap_coefficients(&self) -> (&Field2D<T>, &Field2D<T>, &Field2D<T>, &Field2D<T>) {
        (
            &self.coeffs_u.ap,
            &self.coeffs_u.ap_consistent,
            &self.coeffs_v.ap,
            &self.coeffs_v.ap_consistent,
        )
    }

    fn compute_coefficients_into(
        &mut self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        match component {
            MomentumComponent::U => {
                self.coeffs_u.compute_into(
                    self.grid.nx,
                    self.grid.ny,
                    self.grid.dx,
                    self.grid.dy,
                    dt,
                    component,
                    fields,
                    self.convection_scheme,
                )
            }
            MomentumComponent::V => {
                self.coeffs_v.compute_into(
                    self.grid.nx,
                    self.grid.ny,
                    self.grid.dx,
                    self.grid.dy,
                    dt,
                    component,
                    fields,
                    self.convection_scheme,
                )
            }
        }
    }

    /// Check if a node is on a boundary with Dirichlet BC or Wall BC (which requires Dirichlet treatment)
    fn is_dirichlet_boundary(&self, i: usize, j: usize) -> bool {
        // Check if node is on any boundary with Dirichlet BC or Wall BC (NoSlip/Moving walls)
        if i == 0
            && self.boundary_conditions.get("west").is_some_and(|bc| {
                matches!(
                    bc,
                    BoundaryCondition::Dirichlet { .. } | BoundaryCondition::Wall { .. }
                )
            })
        {
            return true;
        }
        if i == self.grid.nx - 1
            && self.boundary_conditions.get("east").is_some_and(|bc| {
                matches!(
                    bc,
                    BoundaryCondition::Dirichlet { .. } | BoundaryCondition::Wall { .. }
                )
            })
        {
            return true;
        }
        if j == 0
            && self.boundary_conditions.get("south").is_some_and(|bc| {
                matches!(
                    bc,
                    BoundaryCondition::Dirichlet { .. } | BoundaryCondition::Wall { .. }
                )
            })
        {
            return true;
        }
        if j == self.grid.ny - 1
            && self.boundary_conditions.get("north").is_some_and(|bc| {
                matches!(
                    bc,
                    BoundaryCondition::Dirichlet { .. } | BoundaryCondition::Wall { .. }
                )
            })
        {
            return true;
        }
        false
    }

    fn assemble_system_into(
        &mut self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) -> cfd_core::error::Result<()> {
        let coeffs = match component {
            MomentumComponent::U => &self.coeffs_u,
            MomentumComponent::V => &self.coeffs_v,
        };
        let n = self.grid.nx * self.grid.ny;
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let should_rebuild = match component {
            MomentumComponent::U => self.matrix_u.is_none(),
            MomentumComponent::V => self.matrix_v.is_none(),
        };

        if match component {
            MomentumComponent::U => self.rhs_u.as_ref().is_none_or(|r| r.len() != n),
            MomentumComponent::V => self.rhs_v.as_ref().is_none_or(|r| r.len() != n),
        } {
            match component {
                MomentumComponent::U => self.rhs_u = Some(DVector::zeros(n)),
                MomentumComponent::V => self.rhs_v = Some(DVector::zeros(n)),
            }
        } else {
            match component {
                MomentumComponent::U => self.rhs_u.as_mut().unwrap().fill(T::zero()),
                MomentumComponent::V => self.rhs_v.as_mut().unwrap().fill(T::zero()),
            }
        }

        if should_rebuild {
            let mut builder = match component {
                MomentumComponent::U => self.matrix_builder_u.take().unwrap_or_else(|| SparseMatrixBuilder::new(n, n)),
                MomentumComponent::V => self.matrix_builder_v.take().unwrap_or_else(|| SparseMatrixBuilder::new(n, n)),
            };
            let mut rhs = DVector::zeros(n);
            let mask_data = fields.mask.data.as_slice();
            let ap_data = coeffs.ap.data.as_slice();
            let aw_data = coeffs.aw.data.as_slice();
            let ae_data = coeffs.ae.data.as_slice();
            let as_data = coeffs.as_.data.as_slice();
            let an_data = coeffs.an.data.as_slice();
            let source_data = coeffs.source.data.as_slice();
            
            for j in 0..ny {
                for i in 0..nx {
                let idx = j * self.grid.nx + i;

                // Check if this is a masked (solid) cell
                if !mask_data[idx] {
                    // For solid cells: assemble identity equation φ = 0
                    builder.add_entry(idx, idx, T::one())?;
                    // RHS is already zero
                    continue;
                }

                // Check if this is a Dirichlet boundary node
                if self.is_dirichlet_boundary(i, j) {
                    // For Dirichlet BC: assemble identity equation φ = bc_value
                    // This is handled in apply_momentum_boundaries, so skip coefficient assembly
                    builder.add_entry(idx, idx, T::one())?;
                    // RHS will be set by boundary condition handler
                    continue;
                }

                // Central coefficient
                builder.add_entry(idx, idx, ap_data[idx])?;

                // Neighbor coefficients
                if i > 0 {
                    builder.add_entry(idx, idx - 1, -aw_data[idx])?;
                }
                if i < self.grid.nx - 1 {
                    builder.add_entry(idx, idx + 1, -ae_data[idx])?;
                }
                if j > 0 {
                    builder.add_entry(idx, idx - self.grid.nx, -as_data[idx])?;
                }
                if j < self.grid.ny - 1 {
                    builder.add_entry(idx, idx + self.grid.nx, -an_data[idx])?;
                }

                // Source term including pressure gradient
                rhs[idx] = source_data[idx];
            }
        }

        // Apply higher-order boundary conditions for improved near-wall accuracy
        // This provides better velocity gradients near walls for production accuracy
        super::boundary::apply_higher_order_wall_boundaries(
            &mut builder,
            &mut rhs,
            component,
            &self.boundary_conditions,
            &self.grid,
        )?;

        // Apply standard boundary conditions (sets RHS values for Dirichlet, modifies equations for Neumann)
        super::boundary::apply_momentum_boundaries(
            &mut builder,
            &mut rhs,
            component,
            &self.boundary_conditions,
            &self.grid,
        )?;

            match component {
                MomentumComponent::U => {
                    self.matrix_u = Some(builder.build()?);
                    self.matrix_builder_u = Some(SparseMatrixBuilder::new(n, n));
                    self.rhs_u.as_mut().unwrap().copy_from(&rhs);
                }
                MomentumComponent::V => {
                    self.matrix_v = Some(builder.build()?);
                    self.matrix_builder_v = Some(SparseMatrixBuilder::new(n, n));
                    self.rhs_v.as_mut().unwrap().copy_from(&rhs);
                }
            }
        } else {
            let mut matrix = match component {
                MomentumComponent::U => self.matrix_u.take().unwrap(),
                MomentumComponent::V => self.matrix_v.take().unwrap(),
            };
            matrix.values_mut().fill(T::zero());
            
            let mut rhs = DVector::zeros(n);

            let (row_offsets, col_indices, values) = matrix.csr_data_mut();
            let mut update_entry = |r: usize, c: usize, v: T| {
                let start = row_offsets[r];
                let end = row_offsets[r + 1];
                for idx in start..end {
                    if col_indices[idx] == c {
                        values[idx] += v;
                        break;
                    }
                }
            };

            let mask_data = fields.mask.data.as_slice();
            let ap_data = coeffs.ap.data.as_slice();
            let aw_data = coeffs.aw.data.as_slice();
            let ae_data = coeffs.ae.data.as_slice();
            let as_data = coeffs.as_.data.as_slice();
            let an_data = coeffs.an.data.as_slice();
            let source_data = coeffs.source.data.as_slice();

            for j in 0..ny {
                for i in 0..nx {
                    let idx = j * nx + i;

                    if !mask_data[idx] {
                        update_entry(idx, idx, T::one());
                        continue;
                    }

                    if self.is_dirichlet_boundary(i, j) {
                        update_entry(idx, idx, T::one());
                        continue;
                    }

                    update_entry(idx, idx, ap_data[idx]);
                    if i > 0 {
                        update_entry(idx, idx - 1, -aw_data[idx]);
                    }
                    if i < nx - 1 {
                        update_entry(idx, idx + 1, -ae_data[idx]);
                    }
                    if j > 0 {
                        update_entry(idx, idx - nx, -as_data[idx]);
                    }
                    if j < ny - 1 {
                        update_entry(idx, idx + nx, -an_data[idx]);
                    }

                    rhs[idx] = source_data[idx];
                }
            }

            super::boundary::apply_higher_order_wall_boundaries(
                &mut matrix,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;

            super::boundary::apply_momentum_boundaries(
                &mut matrix,
                &mut rhs,
                component,
                &self.boundary_conditions,
                &self.grid,
            )?;
            
            match component {
                MomentumComponent::U => {
                    self.matrix_u = Some(matrix);
                    self.rhs_u.as_mut().unwrap().copy_from(&rhs);
                }
                MomentumComponent::V => {
                    self.matrix_v = Some(matrix);
                    self.rhs_v.as_mut().unwrap().copy_from(&rhs);
                }
            }
        }

        Ok(())
    }

    fn update_velocity(
        &self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        solution: &DVector<T>,
    ) {
        let alpha = self.velocity_relaxation;
        let one_minus_alpha = T::one() - alpha;

        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let idx = j * self.grid.nx + i;
                let computed_value = solution[idx];

                match component {
                    MomentumComponent::U => {
                        if let Some(u) = fields.u.at_mut(i, j) {
                            let old_value = *u;
                            // Enforce zero velocity in masked (solid) cells
                            if fields.mask.at(i, j) {
                                // Under-relaxation: u_new = α * u_computed + (1-α) * u_old
                                *u = alpha * computed_value + one_minus_alpha * old_value;
                            } else {
                                *u = T::zero();
                            }
                        }
                    }
                    MomentumComponent::V => {
                        if let Some(v) = fields.v.at_mut(i, j) {
                            let old_value = *v;
                            // Enforce zero velocity in masked (solid) cells
                            if fields.mask.at(i, j) {
                                // Under-relaxation: v_new = α * v_computed + (1-α) * v_old
                                *v = alpha * computed_value + one_minus_alpha * old_value;
                            } else {
                                *v = T::zero();
                            }
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::StructuredGrid2D;
    use cfd_core::physics::boundary::{BoundaryCondition, WallType};
    use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
    use cfd_math::linear_solver::IterativeLinearSolver;
    use nalgebra::Vector3;

    #[test]
    fn moving_north_wall_assembles_nonzero_coupling_to_interior() {
        let grid = StructuredGrid2D::new(4, 4, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)
            .expect("grid creation failed");
        let mut solver = MomentumSolver::new(&grid);
        solver.set_boundary(
            "north".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::Moving {
                    velocity: Vector3::new(1.0, 0.0, 0.0),
                },
            },
        );
        solver.set_boundary(
            "south".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );
        solver.set_boundary(
            "east".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );
        solver.set_boundary(
            "west".to_string(),
            BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );

        let mut fields = SimulationFields::new(4, 4);
        fields.density.map_inplace(|d| *d = 1.0);
        fields.viscosity.map_inplace(|v| *v = 0.01);

        solver
            .solve_with_coefficients(MomentumComponent::U, &mut fields, 0.0125)
            .expect("momentum solve failed");

        let matrix = solver
            .matrix_u
            .as_ref()
            .expect("momentum matrix should be assembled");
        let rhs = solver.rhs_u.as_ref().expect("momentum rhs should be assembled");
        let row = 2 * grid.nx + 2;
        let start = matrix.row_offsets()[row];
        let end = matrix.row_offsets()[row + 1];
        let cols = &matrix.col_indices()[start..end];
        let vals = &matrix.values()[start..end];

        let top_row = (grid.ny - 1) * grid.nx + 2;
        let top_coeff = cols
            .iter()
            .zip(vals.iter())
            .find(|(col, _)| **col == top_row)
            .map(|(_, val)| *val);

        assert!(
            top_coeff.is_some(),
            "interior momentum row must couple to the moving lid row"
        );
        assert!(
            top_coeff.expect("checked above").abs() > 0.0,
            "moving lid coupling coefficient must be non-zero"
        );

        let direct_solution = cfd_math::linear_solver::DirectSparseSolver::default()
            .solve(matrix, rhs)
            .expect("direct sparse solve should succeed");
        let mut gmres_solution = DVector::zeros(matrix.nrows());
        solver
            .linear_solver
            .solve(
                matrix,
                rhs,
                &mut gmres_solution,
                None::<&IdentityPreconditioner>,
            )
            .expect("GMRES solve should succeed");
        assert!(
            direct_solution[row] > 0.0,
            "direct solve should produce a positive interior response"
        );
        assert!(
            gmres_solution[row] > 0.0,
            "GMRES solve should also produce a positive interior response"
        );
    }
}
