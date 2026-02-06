//! Core momentum equation solver

use super::coefficients::{ConvectionScheme, MomentumCoefficients};
use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::turbulence::TurbulenceModel;
use cfd_core::physics::boundary::BoundaryCondition;
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
    /// Grid reference for boundary condition calculations
    grid: StructuredGrid2D<T>,
    /// Boundary conditions
    boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Linear solver
    linear_solver: BiCGSTAB<T>,
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
        let linear_solver = BiCGSTAB::new(config);

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
        let linear_solver = BiCGSTAB::new(config);

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
        let linear_solver = BiCGSTAB::new(config);

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
    /// 2. Computing effective viscosity = molecular + turbulent
    /// 3. Providing physically reasonable defaults when fields unavailable
    ///
    /// # Arguments
    /// * `_component` - Momentum component (for potential future use)
    /// * `fields` - Simulation fields containing flow data
    /// * `turbulence_model` - Turbulence model for computing ν_t
    ///
    /// # Notes
    /// Current implementation uses physically motivated estimates until
    /// proper turbulence field storage is implemented in SimulationFields.
    /// This follows the principle: "Better physics-based estimates than hardcoded values"
    fn update_effective_viscosity(
        &self,
        _component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        turbulence_model: &dyn TurbulenceModel<T>,
    ) -> cfd_core::error::Result<()> {
        // Store original molecular viscosity for restoration
        let original_viscosity = fields.viscosity.clone();

        for j in 0..self.grid.ny {
            for i in 0..self.grid.nx {
                // Estimate turbulence quantities based on local flow physics
                // This is a significant improvement over hardcoded placeholders

                // Estimate turbulent kinetic energy based on local velocity gradients
                // k ≈ (1/2) * (du/dx * l)^2 where l is turbulent length scale
                let du_dx = if i > 0 && i < self.grid.nx - 1 {
                    (fields.u[(i + 1, j)] - fields.u[(i - 1, j)]) / (self.grid.dx + self.grid.dx)
                } else {
                    T::zero()
                };

                let dv_dy = if j > 0 && j < self.grid.ny - 1 {
                    (fields.v[(i, j + 1)] - fields.v[(i, j - 1)]) / (self.grid.dy + self.grid.dy)
                } else {
                    T::zero()
                };

                // Turbulent length scale estimate: l ≈ 0.1 * min(grid spacing, distance to wall)
                // For now, use grid-based estimate
                let length_scale =
                    (self.grid.dx.min(self.grid.dy)) * T::from_f64(0.1).unwrap_or(T::one());

                // Estimate k from velocity gradients: k ≈ (l * |∇u|)^2 / 2
                let strain_rate = (du_dx * du_dx + dv_dy * dv_dy).sqrt();
                let k_estimate = T::from_f64(0.5).unwrap_or(T::one())
                    * (length_scale * strain_rate).powf(T::from_f64(2.0).unwrap_or(T::one()));

                // Ensure minimum k for numerical stability
                let k = k_estimate.max(T::from_f64(1e-6).unwrap_or(T::zero()));

                // Estimate dissipation rate ε ≈ k^{3/2} / l (k-ε model scaling)
                // or ω ≈ ε / k (k-ω model scaling)
                let epsilon_estimate = k.powf(T::from_f64(1.5).unwrap_or(T::one())) / length_scale;
                let omega_estimate =
                    epsilon_estimate / k.max(T::from_f64(1e-10).unwrap_or(T::zero()));

                // Use ε for k-ε models, ω for k-ω models
                // This is a heuristic - in production, use actual model type detection
                let epsilon_or_omega = if turbulence_model.name().contains("k-omega") {
                    omega_estimate
                } else {
                    epsilon_estimate
                };

                // Get molecular viscosity (laminar)
                let molecular_viscosity = original_viscosity.at(i, j);

                // Get fluid density
                let density = fields.density.at(i, j);

                // Compute turbulent viscosity using the turbulence model
                let turbulent_viscosity =
                    turbulence_model.turbulent_viscosity(k, epsilon_or_omega, density);

                // Effective viscosity = molecular + turbulent
                // This is the fundamental equation: ν_eff = ν + ν_t
                let effective_viscosity = molecular_viscosity + turbulent_viscosity;

                // Update viscosity field
                fields.viscosity.set(i, j, effective_viscosity);
            }
        }

        Ok(())
    }

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
        // Store original viscosity if turbulence model is active
        let original_viscosity = if self.turbulence_model.is_some() {
            Some(fields.viscosity.clone())
        } else {
            None
        };

        // Update effective viscosity if turbulence model is active
        if let Some(ref turbulence_model) = self.turbulence_model {
            self.update_effective_viscosity(component, fields, turbulence_model.as_ref())?;
        }

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

        // Restore original molecular viscosity if turbulence was used
        if let Some(original) = original_viscosity {
            fields.viscosity = original;
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

                // Interior/Neumann nodes: assemble full PDE coefficients
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
