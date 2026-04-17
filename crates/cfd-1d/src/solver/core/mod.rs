//! Modularized network solver for 1D CFD analysis
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies.

mod anderson_acceleration;
mod convergence;
mod geometry;
mod linear_system;
mod matrix_assembly;
mod problem;
mod solver_detection;
mod state;
mod status;
/// Transient solvers for composition and droplet tracking over time.
pub mod transient;
/// Workspace state and allocation for the 1D solver loop.
pub mod workspace;

pub use convergence::ConvergenceChecker;
pub use geometry::NetworkDomain;
pub use linear_system::{LinearSolverMethod, LinearSystemSolver};
pub use matrix_assembly::MatrixAssembler;
pub use problem::NetworkProblem;
pub use state::NetworkState;
pub use status::{PrimarySolveDiagnostics, PrimarySolveError, SolveFailureReason, SolvePathStatus};
pub use workspace::SolverWorkspace;

pub use transient::composition::{
    BloodEdgeTransportConfig, CompositionState, EdgeFlowEvent, InletCompositionEvent,
    InletHematocritEvent, MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig,
    TransientCompositionSimulator, BLOOD_PLASMA_FLUID_ID, BLOOD_RBC_FLUID_ID,
};
pub use transient::droplets::{
    ChannelOccupancy, DropletBoundary, DropletInjection, DropletPosition, DropletSnapshot,
    DropletSplitPolicy, DropletState, DropletTrackingState, SplitMode, TransientDropletSimulator,
};

use crate::domain::network::Network;
use cfd_core::compute::solver::{Configurable, Solver, Validatable};
use cfd_core::error::Result;
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind};
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// Solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig<T: RealField + Copy> {
    /// Convergence tolerance for solution accuracy
    pub tolerance: T,
    /// Maximum number of solver iterations before termination
    pub max_iterations: usize,
}

impl<T: RealField + Copy> cfd_core::compute::solver::SolverConfiguration<T> for SolverConfig<T> {
    fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }

    fn use_preconditioning(&self) -> bool {
        false // No preconditioning for network solver
    }
}

/// Main network solver implementing the core CFD suite trait system
pub struct NetworkSolver<T: RealField + Copy, F: FluidTrait<T> = ConstantPropertyFluid<T>> {
    /// Solver configuration
    config: SolverConfig<T>,
    /// Matrix assembler for building the linear system
    assembler: MatrixAssembler<T>,
    /// Convergence checker
    convergence: ConvergenceChecker<T>,
    _phantom: std::marker::PhantomData<F>,
}

impl<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone> Default
    for NetworkSolver<T, F>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone> NetworkSolver<T, F> {
    /// Create a new network solver with default configuration
    #[must_use]
    pub fn new() -> Self {
        // Default tolerance: 1e-6, falling back to a multiple of machine epsilon
        // so that the solver never panics on exotic numeric types.
        let tolerance = T::from_f64(1e-6).expect("Mathematical constant conversion compromised");
        let config = SolverConfig {
            tolerance,
            max_iterations: 1000,
        };
        let convergence = ConvergenceChecker::new(config.tolerance);
        Self {
            config,
            assembler: MatrixAssembler::new(),
            convergence,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create with specific configuration
    pub fn with_config(config: SolverConfig<T>) -> Self {
        Self {
            assembler: MatrixAssembler::new(),
            convergence: ConvergenceChecker::new(config.tolerance),
            config,
            _phantom: std::marker::PhantomData,
        }
    }

    fn validate_network_contract(&self, network: &Network<T, F>) -> Result<()> {
        if network.node_count() == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has no nodes".to_string(),
            ));
        }
        if self.config.tolerance <= T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Tolerance must be positive".to_string(),
            ));
        }
        network.validate_coefficients()?;
        for props in network.properties.values() {
            if props.length <= T::zero() || !props.length.is_finite() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Edge '{}' has invalid physical length",
                    props.id
                )));
            }
            if props.area <= T::zero() || !props.area.is_finite() {
                return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                    "Edge '{}' has invalid cross-sectional area",
                    props.id
                )));
            }
            if let Some(d_h) = props.hydraulic_diameter {
                if d_h <= T::zero() || !d_h.is_finite() {
                    return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                        "Edge '{}' has invalid hydraulic diameter",
                        props.id
                    )));
                }
            }
        }
        Ok(())
    }
}

impl<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    > NetworkSolver<T, F>
{
    /// Solve the network flow problem iteratively for non-linear systems.
    ///
    /// ## Algorithm: Anderson-Accelerated Picard Iteration
    ///
    /// **Theorem** (Walker & Ni 2011): Anderson acceleration applied to a
    /// fixed-point iteration x_{k+1} = G(x_k) achieves superlinear convergence
    /// for contractive mappings, reducing the Picard iteration count by a factor
    /// of 2–5× for typical microfluidic networks.
    ///
    /// The conductance-matrix Laplacian from `MatrixAssembler` is structurally
    /// invariant across Picard iterations — only the edge conductance magnitudes
    /// change as flow-dependent friction factors are updated. This means:
    ///
    /// 1. **SPD detection is hoisted to the first iteration**: the sign structure
    ///    of the Laplacian (positive diagonal, non-positive off-diagonal) is
    ///    preserved across conductance updates, so the solver method selection
    ///    (CG vs BiCGSTAB) is performed once.
    ///
    /// 2. **Anderson acceleration with depth m=5** is applied to the pressure
    ///    solution vector, treating each Picard iterate as a fixed-point map.
    ///    The MGS-QR subproblem (from `cfd-math`) provides O(m²n) per-step
    ///    cost with guaranteed numerical stability.
    ///
    /// 3. **The linear solver object is allocated once** and reused across
    ///    iterations, eliminating per-iteration `LinearSystemSolver::new()`.
    pub fn solve_network(&self, problem: &NetworkProblem<T, F>) -> Result<Network<T, F>> {
        self.solve_network_with_diagnostics(problem)
            .map(|(network, _)| network)
            .map_err(PrimarySolveError::into_source)
    }

    /// Solve a network directly, consuming ownership to avoid cloning the graph.
    ///
    /// This path is intended for transient workflows that repeatedly mutate and
    /// re-solve the same network state. It preserves the same validation and
    /// diagnostics behavior as [`Self::solve_network`] while eliminating the
    /// internal clone required by the borrowed `NetworkProblem` path.
    pub fn solve_owned_network(&self, network: Network<T, F>) -> Result<Network<T, F>> {
        self.solve_owned_network_with_diagnostics(network)
            .map(|(network, _)| network)
            .map_err(PrimarySolveError::into_source)
    }

    /// Solve the network and return explicit diagnostics for the trusted primary path.
    #[allow(clippy::too_many_lines, clippy::result_large_err)]
    pub fn solve_network_with_diagnostics(
        &self,
        problem: &NetworkProblem<T, F>,
    ) -> std::result::Result<(Network<T, F>, PrimarySolveDiagnostics), PrimarySolveError> {
        self.solve_owned_network_with_diagnostics(problem.network.clone())
    }

    /// Solve the network and return explicit diagnostics without cloning the input network.
    #[allow(clippy::too_many_lines, clippy::result_large_err)]
    pub fn solve_owned_network_with_diagnostics(
        &self,
        network: Network<T, F>,
    ) -> std::result::Result<(Network<T, F>, PrimarySolveDiagnostics), PrimarySolveError> {
        self.validate_network_contract(&network).map_err(|source| {
            PrimarySolveError::new(
                SolveFailureReason::InvalidGeometryContract,
                PrimarySolveDiagnostics::default(),
                source,
            )
        })?;

        let mut network = network;
        let n = network.node_count();
        let anderson_depth = 5;
        let mut workspace = workspace::SolverWorkspace::new(n, anderson_depth);

        MatrixAssembler::<T>::classify_boundary_conditions_into(
            &network,
            &mut workspace.dirichlet_values,
            &mut workspace.neumann_sources,
        )
        .map_err(|source| {
            PrimarySolveError::new(
                SolveFailureReason::InvalidGeometryContract,
                PrimarySolveDiagnostics::default(),
                source,
            )
        })?;

        let mut last_flow_rates = network.flow_rates.clone();
        let mut diagnostics = PrimarySolveDiagnostics::default();
        network.residuals.clear();
        network.residuals.reserve(self.config.max_iterations);

        if Self::is_linear_static_network(&network) {
            diagnostics.matrix_treated_as_linear_static = true;
            let matrix = self
                .assembler
                .assemble_into(&network, &mut workspace)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::MatrixAssemblyInvalid,
                        diagnostics.clone(),
                        source,
                    )
                })?;
            Self::validate_linear_system(&matrix, &workspace.rhs).map_err(|source| {
                PrimarySolveError::new(
                    SolveFailureReason::MatrixAssemblyInvalid,
                    diagnostics.clone(),
                    source,
                )
            })?;
            let method = Self::detect_solver_method(&matrix);
            diagnostics.linear_solver_method = Some(method);
            diagnostics.picard_iterations = 1;
            network.last_solver_method = Some(method);
            workspace.linear_solution.fill(T::zero());
            let solution = LinearSystemSolver::new()
                .with_method(method)
                .with_tolerance(self.config.tolerance)
                .with_max_iterations(self.config.max_iterations)
                .solve_with_initial_guess(&matrix, &workspace.rhs, &mut workspace.linear_solution)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::LinearSolverFailure,
                        diagnostics.clone(),
                        source,
                    )
                })?;
            if !Self::vector_is_finite(&solution) {
                return Err(PrimarySolveError::new(
                    SolveFailureReason::NonFiniteResidual,
                    diagnostics,
                    Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: "non-finite linear-static solution".to_string(),
                    }),
                ));
            }
            let residual_norm = Self::compute_residual_norm(&matrix, &solution, &workspace.rhs, n);
            diagnostics.last_residual_norm = Self::scalar_to_f64(residual_norm);
            diagnostics.last_solution_change_norm = Self::scalar_to_f64(solution.norm());
            if !residual_norm.is_finite() {
                return Err(PrimarySolveError::new(
                    SolveFailureReason::NonFiniteResidual,
                    diagnostics,
                    Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: "non-finite linear-static residual".to_string(),
                    }),
                ));
            }
            self.update_network_solution(&mut network, &solution)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::MatrixAssemblyInvalid,
                        diagnostics.clone(),
                        source,
                    )
                })?;
            return Ok((network, diagnostics));
        }

        let mut selected_method: Option<LinearSolverMethod> = None;
        let mut adaptive_solver: Option<LinearSystemSolver<T>> = None;

        for iter in 0..self.config.max_iterations {
            diagnostics.picard_iterations = iter + 1;
            let matrix = self
                .assembler
                .assemble_into(&network, &mut workspace)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::MatrixAssemblyInvalid,
                        diagnostics.clone(),
                        source,
                    )
                })?;
            Self::validate_linear_system(&matrix, &workspace.rhs).map_err(|source| {
                PrimarySolveError::new(
                    SolveFailureReason::MatrixAssemblyInvalid,
                    diagnostics.clone(),
                    source,
                )
            })?;

            if selected_method.is_none() {
                let method = Self::detect_solver_method(&matrix);
                selected_method = Some(method);
                network.last_solver_method = Some(method);
                diagnostics.linear_solver_method = Some(method);
                adaptive_solver = Some(
                    LinearSystemSolver::new()
                        .with_method(method)
                        .with_tolerance(self.config.tolerance)
                        .with_max_iterations(self.config.max_iterations),
                );
            }

            workspace
                .linear_solution
                .copy_from(&workspace.last_solution);
            let picard_solution = adaptive_solver
                .as_ref()
                .expect("adaptive solver must be initialized after method detection")
                .solve_with_initial_guess(&matrix, &workspace.rhs, &mut workspace.linear_solution)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::LinearSolverFailure,
                        diagnostics.clone(),
                        source,
                    )
                })?;
            if !Self::vector_is_finite(&picard_solution) {
                return Err(PrimarySolveError::new(
                    SolveFailureReason::NonFiniteResidual,
                    diagnostics,
                    Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: "non-finite Picard iterate".to_string(),
                    }),
                ));
            }
            let picard_residual =
                Self::compute_residual_norm(&matrix, &picard_solution, &workspace.rhs, n);
            if !picard_residual.is_finite() {
                diagnostics.last_residual_norm = Self::scalar_to_f64(picard_residual);
                return Err(PrimarySolveError::new(
                    SolveFailureReason::NonFiniteResidual,
                    diagnostics,
                    Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: "non-finite Picard residual".to_string(),
                    }),
                ));
            }

            let solution = Self::select_next_iterate(
                iter,
                n,
                picard_solution,
                &mut workspace,
                anderson_depth,
                &matrix,
                picard_residual,
            );

            let residual_norm = Self::compute_residual_norm(&matrix, &solution, &workspace.rhs, n);
            let rhs_norm = workspace.rhs.norm();
            let solution_change_norm = (&solution - &workspace.last_solution).norm();
            diagnostics.last_residual_norm = Self::scalar_to_f64(residual_norm);
            diagnostics.last_solution_change_norm = Self::scalar_to_f64(solution_change_norm);
            if !residual_norm.is_finite() || !solution_change_norm.is_finite() {
                return Err(PrimarySolveError::new(
                    SolveFailureReason::NonFiniteResidual,
                    diagnostics,
                    Error::Numerical(NumericalErrorKind::InvalidValue {
                        value: "non-finite primary residual or solution delta".to_string(),
                    }),
                ));
            }
            network.residuals.push(residual_norm);

            let converged = self
                .convergence
                .has_converged_dual(&solution, &workspace.last_solution, residual_norm, rhs_norm)
                .map_err(|source| {
                    let reason = match &source {
                        Error::Convergence(ConvergenceErrorKind::InvalidValue)
                        | Error::Numerical(NumericalErrorKind::InvalidValue { .. }) => {
                            SolveFailureReason::NonFiniteResidual
                        }
                        _ => SolveFailureReason::MaxIterationsExceeded,
                    };
                    PrimarySolveError::new(reason, diagnostics.clone(), source)
                })?;

            self.update_network_solution(&mut network, &solution)
                .map_err(|source| {
                    PrimarySolveError::new(
                        SolveFailureReason::MatrixAssemblyInvalid,
                        diagnostics.clone(),
                        source,
                    )
                })?;

            let mut flow_diff_sq = T::zero();
            let mut flow_norm_sq = T::zero();
            for (i, &new_flow) in network.flow_rates.iter().enumerate() {
                let old_flow = last_flow_rates.get(i).copied().unwrap_or(T::zero());
                let diff = new_flow - old_flow;
                flow_diff_sq += diff * diff;
                flow_norm_sq += new_flow * new_flow;
            }
            let flow_change = <T as Float>::sqrt(flow_diff_sq);
            let flow_norm = <T as Float>::sqrt(flow_norm_sq);
            let relative_flow_change = if flow_norm > T::default_epsilon() {
                flow_change / flow_norm
            } else {
                flow_change
            };
            let flows_converged = relative_flow_change < self.config.tolerance;

            if converged && flows_converged {
                return Ok((network, diagnostics));
            }

            workspace.last_solution = solution;
            last_flow_rates.clone_from(&network.flow_rates);
        }

        Err(PrimarySolveError::new(
            SolveFailureReason::MaxIterationsExceeded,
            diagnostics,
            Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            }),
        ))
    }

    fn update_network_solution(
        &self,
        network: &mut Network<T, F>,
        solution: &nalgebra::DVector<T>,
    ) -> Result<()> {
        network.update_from_solution(solution)
    }
}

impl<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    > Solver<T> for NetworkSolver<T, F>
{
    type Problem = NetworkProblem<T, F>;
    type Solution = Network<T, F>;

    fn solve(&self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_network(problem)
    }

    fn name(&self) -> &'static str {
        "NetworkSolver"
    }
}

impl<T: RealField + Copy, F: FluidTrait<T>> Configurable<T> for NetworkSolver<T, F> {
    type Config = SolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn set_config(&mut self, config: Self::Config) {
        self.config = config;
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive, F: FluidTrait<T> + Clone> Validatable<T>
    for NetworkSolver<T, F>
{
    type Problem = NetworkProblem<T, F>;

    fn validate_problem(&self, problem: &Self::Problem) -> Result<()> {
        self.validate_network_contract(&problem.network)
    }
}
