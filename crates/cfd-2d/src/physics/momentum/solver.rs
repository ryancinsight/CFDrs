//! Core momentum equation solver

use super::coefficients::{ConvectionScheme, MomentumCoefficients};
use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::turbulence::TurbulenceModel;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::IterativeSolverConfig;
use cfd_math::linear_solver::{GMRES, IterativeLinearSolver};

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
    /// Last computed A_p coefficients for U
    ap_u: Field2D<T>,
    /// Last computed A_p_consistent coefficients for U (SIMPLEC)
    ap_consistent_u: Field2D<T>,
    /// Last computed A_p coefficients for V
    ap_v: Field2D<T>,
    /// Last computed A_p_consistent coefficients for V (SIMPLEC)
    ap_consistent_v: Field2D<T>,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::ToPrimitive> MomentumSolver<T> {
    /// Create new momentum solver with default deferred correction scheme
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        // Use relaxed tolerance for momentum equations in SIMPLE algorithms
        // Momentum equations don't need extreme accuracy - they get corrected by pressure
        let config = IterativeSolverConfig {
            max_iterations: 2000, // Reduced for efficiency
            tolerance: T::from_f64(5e-2).unwrap_or_else(|| T::from_f64(1e-3).unwrap()), // Very relaxed tolerance for SIMPLE
            use_preconditioner: true,
            use_parallel_spmv: false,
        };
        let linear_solver = GMRES::new(config, 30);



        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(T::one),
            turbulence_model: None,
            ap_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_v: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_v: Field2D::new(grid.nx, grid.ny, T::one()),
        }
    }

    /// Create new momentum solver with specified convection scheme
    pub fn with_convection_scheme(grid: &StructuredGrid2D<T>, scheme: ConvectionScheme) -> Self {
        let config = IterativeSolverConfig::default();
        let linear_solver = GMRES::new(config, 30);



        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: scheme,
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(T::one),
            turbulence_model: None,
            ap_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_v: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_v: Field2D::new(grid.nx, grid.ny, T::one()),
        }
    }

    /// Create new momentum solver with parallel SpMV enabled for multi-core performance
    #[must_use]
    pub fn with_parallel_spmv(grid: &StructuredGrid2D<T>) -> Self {
        let config = IterativeSolverConfig::default().with_parallel_spmv();
        let linear_solver = GMRES::new(config, 30);



        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(T::one),
            turbulence_model: None,
            ap_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_u: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_v: Field2D::new(grid.nx, grid.ny, T::one()),
            ap_consistent_v: Field2D::new(grid.nx, grid.ny, T::one()),
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

    /// Update effective viscosity field using turbulence model
    ///
    /// This method properly integrates with the turbulence model by:
    /// 1. Using estimated turbulence quantities based on flow physics

    /// Solve momentum equation for specified component
    pub fn solve(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        let _ = self.solve_with_coefficients(component, fields, dt)?;
        Ok(())
    }

    /// Solve momentum equation and return coefficients for Rhie-Chow interpolation
    pub fn solve_with_coefficients(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<MomentumCoefficients<T>> {


        // Compute coefficients
        let coeffs = self.compute_coefficients(component, fields, dt)?;

        // Diagnostic: Check if coefficients are non-zero (helps debug false convergence)
        #[cfg(debug_assertions)]
        {
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
        let rhs_norm = rhs.norm();
        self.linear_solver.solve(
            &matrix,
            &rhs,
            &mut solution,
            None::<&IdentityPreconditioner>,
        )?;



        // Update velocity field
        self.update_velocity(component, fields, &solution);

        // Store A_p and A_p_consistent coefficients for pressure correction
        match component {
            MomentumComponent::U => {
                self.ap_u = coeffs.ap.clone();
                self.ap_consistent_u = coeffs.ap_consistent.clone();
            }
            MomentumComponent::V => {
                self.ap_v = coeffs.ap.clone();
                self.ap_consistent_v = coeffs.ap_consistent.clone();
            }
        }



        Ok(coeffs)
    }

    /// Get the last computed A_p and A_p_consistent coefficients for both velocity components
    pub fn get_ap_coefficients(&self) -> (Field2D<T>, Field2D<T>, Field2D<T>, Field2D<T>) {
        (
            self.ap_u.clone(),
            self.ap_consistent_u.clone(),
            self.ap_v.clone(),
            self.ap_consistent_v.clone(),
        )
    }

    fn compute_coefficients(
        &self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<MomentumCoefficients<T>> {
        MomentumCoefficients::compute(
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

    fn assemble_system(
        &self,
        coeffs: &MomentumCoefficients<T>,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) -> cfd_core::error::Result<(SparseMatrix<T>, DVector<T>)> {
        let n = self.grid.nx * self.grid.ny;
        let mut builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::zeros(n);

        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                let idx = j * self.grid.nx + i;

                // Check if this is a masked (solid) cell
                if !fields.mask.at(i, j) {
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
                builder.add_entry(idx, idx, coeffs.ap.at(i, j))?;

                // Neighbor coefficients
                if i > 0 {
                    builder.add_entry(idx, idx - 1, -coeffs.aw.at(i, j))?;
                }
                if i < self.grid.nx - 1 {
                    builder.add_entry(idx, idx + 1, -coeffs.ae.at(i, j))?;
                }
                if j > 0 {
                    builder.add_entry(idx, idx - self.grid.nx, -coeffs.as_.at(i, j))?;
                }
                if j < self.grid.ny - 1 {
                    builder.add_entry(idx, idx + self.grid.nx, -coeffs.an.at(i, j))?;
                }

                // Source term including pressure gradient
                rhs[idx] = coeffs.source.at(i, j);
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
