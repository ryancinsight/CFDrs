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
use cfd_math::linear_solver::{IterativeSolverConfig, GMRES};

use crate::scalar::Cfd2dScalar;
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use eunomia::FloatElement;
use leto::Array1;
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
pub struct MomentumSolver<T: Cfd2dScalar + Copy> {
    /// Grid reference for boundary condition calculations
    pub(crate) grid: StructuredGrid2D<T>,
    /// Boundary conditions
    pub(crate) boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Linear solver
    pub(crate) linear_solver: GMRES<T>,

    /// Convection discretization scheme
    pub(crate) convection_scheme: ConvectionScheme,
    /// Velocity under-relaxation factor (0 < α ≤ 1, default 0.7)
    /// u_new = α * u_computed + (1-α) * u_old
    pub(crate) velocity_relaxation: T,
    /// Optional turbulence model for computing turbulent viscosity
    pub(crate) turbulence_model: Option<Box<dyn TurbulenceModel<T>>>,

    /// U-momentum coefficients
    pub(crate) coeffs_u: MomentumCoefficients<T>,
    /// V-momentum coefficients
    pub(crate) coeffs_v: MomentumCoefficients<T>,

    pub(crate) matrix_u: Option<SparseMatrix<T>>,
    pub(crate) matrix_v: Option<SparseMatrix<T>>,
    pub(crate) matrix_builder_u: Option<SparseMatrixBuilder<T>>,
    pub(crate) matrix_builder_v: Option<SparseMatrixBuilder<T>>,
    pub(crate) rhs_u: Option<Array1<T>>,
    pub(crate) rhs_v: Option<Array1<T>>,
}

impl<T: Cfd2dScalar + Copy + FloatElement> MomentumSolver<T> {
    /// Create new momentum solver with default deferred correction scheme
    pub fn new(grid: &StructuredGrid2D<T>) -> Self {
        let linear_solver = GMRES::new(Self::linear_solver_config(), 30);

        Self {
            grid: grid.clone(),
            boundary_conditions: HashMap::new(),
            linear_solver,
            convection_scheme: ConvectionScheme::default(),
            velocity_relaxation: <T as FloatElement>::from_f64(0.7),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::U,
                &SimulationFields::new(grid.nx, grid.ny),
                ConvectionScheme::default(),
            )
            .unwrap(),
            coeffs_v: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::V,
                &SimulationFields::new(grid.nx, grid.ny),
                ConvectionScheme::default(),
            )
            .unwrap(),
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
            velocity_relaxation: <T as FloatElement>::from_f64(0.7),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::U,
                &SimulationFields::new(grid.nx, grid.ny),
                scheme,
            )
            .unwrap(),
            coeffs_v: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::V,
                &SimulationFields::new(grid.nx, grid.ny),
                scheme,
            )
            .unwrap(),
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
            velocity_relaxation: <T as FloatElement>::from_f64(0.7),
            turbulence_model: None,
            coeffs_u: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::U,
                &SimulationFields::new(grid.nx, grid.ny),
                ConvectionScheme::default(),
            )
            .unwrap(),
            coeffs_v: MomentumCoefficients::compute(
                grid.nx,
                grid.ny,
                T::one(),
                T::one(),
                T::one(),
                MomentumComponent::V,
                &SimulationFields::new(grid.nx, grid.ny),
                ConvectionScheme::default(),
            )
            .unwrap(),
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

    fn linear_solver_config() -> IterativeSolverConfig<T> {
        IterativeSolverConfig {
            max_iterations: crate::constants::solver::DEFAULT_MAX_ITERATIONS,
            tolerance: <T as FloatElement>::from_f64(crate::constants::solver::DEFAULT_TOLERANCE),
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

    pub(crate) fn compute_coefficients_into(
        &mut self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        dt: T,
    ) -> cfd_core::error::Result<()> {
        match component {
            MomentumComponent::U => self.coeffs_u.compute_into(
                self.grid.nx,
                self.grid.ny,
                self.grid.dx,
                self.grid.dy,
                dt,
                component,
                fields,
                self.convection_scheme,
            ),
            MomentumComponent::V => self.coeffs_v.compute_into(
                self.grid.nx,
                self.grid.ny,
                self.grid.dx,
                self.grid.dy,
                dt,
                component,
                fields,
                self.convection_scheme,
            ),
        }
    }

    /// Check if a node is on a boundary with Dirichlet BC or Wall BC (which requires Dirichlet treatment)
    pub(crate) fn is_dirichlet_boundary(
        &self,
        i: usize,
        j: usize,
        component: MomentumComponent,
    ) -> bool {
        // South boundary
        if j == 0 {
            if let Some(bc) = self.boundary_conditions.get("south") {
                return match bc {
                    BoundaryCondition::Dirichlet { .. }
                    | BoundaryCondition::VelocityInlet { .. }
                    | BoundaryCondition::CharacteristicInlet { .. } => true,
                    BoundaryCondition::Wall { wall_type } => match wall_type {
                        cfd_core::physics::boundary::WallType::Slip => {
                            matches!(component, MomentumComponent::V)
                        }
                        _ => true,
                    },
                    BoundaryCondition::Symmetry => matches!(component, MomentumComponent::V),
                    _ => false,
                };
            }
        }
        // North boundary
        if j == self.grid.ny - 1 {
            if let Some(bc) = self.boundary_conditions.get("north") {
                return match bc {
                    BoundaryCondition::Dirichlet { .. }
                    | BoundaryCondition::VelocityInlet { .. }
                    | BoundaryCondition::CharacteristicInlet { .. } => true,
                    BoundaryCondition::Wall { wall_type } => match wall_type {
                        cfd_core::physics::boundary::WallType::Slip => {
                            matches!(component, MomentumComponent::V)
                        }
                        _ => true,
                    },
                    BoundaryCondition::Symmetry => matches!(component, MomentumComponent::V),
                    _ => false,
                };
            }
        }
        // West boundary
        if i == 0 {
            if let Some(bc) = self.boundary_conditions.get("west") {
                return match bc {
                    BoundaryCondition::Dirichlet { .. }
                    | BoundaryCondition::VelocityInlet { .. }
                    | BoundaryCondition::CharacteristicInlet { .. } => true,
                    BoundaryCondition::Wall { wall_type } => match wall_type {
                        cfd_core::physics::boundary::WallType::Slip => {
                            matches!(component, MomentumComponent::U)
                        }
                        _ => true,
                    },
                    BoundaryCondition::Symmetry => matches!(component, MomentumComponent::U),
                    _ => false,
                };
            }
        }
        // East boundary
        if i == self.grid.nx - 1 {
            if let Some(bc) = self.boundary_conditions.get("east") {
                return match bc {
                    BoundaryCondition::Dirichlet { .. }
                    | BoundaryCondition::VelocityInlet { .. }
                    | BoundaryCondition::CharacteristicInlet { .. } => true,
                    BoundaryCondition::Wall { wall_type } => match wall_type {
                        cfd_core::physics::boundary::WallType::Slip => {
                            matches!(component, MomentumComponent::U)
                        }
                        _ => true,
                    },
                    BoundaryCondition::Symmetry => matches!(component, MomentumComponent::U),
                    _ => false,
                };
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::StructuredGrid2D;
    use cfd_core::physics::boundary::{BoundaryCondition, WallType};
    use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
    use cfd_math::linear_solver::{DirectSparseSolver, IterativeLinearSolver};
    use leto::geometry::Vector3;
    use leto::Array1;

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
        let rhs = solver
            .rhs_u
            .as_ref()
            .expect("momentum rhs should be assembled");
        let row = 2 * grid.nx + 2;
        let start = matrix.row_ptr()[row];
        let end = matrix.row_ptr()[row + 1];
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

        let direct_solution = DirectSparseSolver::default()
            .solve(matrix, rhs)
            .expect("direct sparse solve should succeed");
        let mut gmres_solution = Array1::from_elem([matrix.nrows()], 0.0);
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
