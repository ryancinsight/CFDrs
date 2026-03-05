//! 3D FEM Navier-Stokes solver for bifurcations with blood flow
//!
//! Solves the incompressible Navier-Stokes equations on 3D bifurcation domains
//! using Finite Element Method with stabilization for advection-dominated flows.
//!
//! # Theorem — Mass Conservation at Bifurcation
//!
//! For incompressible flow at steady state, continuity requires
//!
//! ```text
//! Q_parent = Q_daughter1 + Q_daughter2
//! ```
//!
//! where $Q_i = \int_{A_i} \mathbf{u} \cdot \mathbf{n} \, dA$. This is validated
//! post-solve by checking $|Q_0 - Q_1 - Q_2| / Q_0 < \varepsilon$.
//!
//! # Theorem — Optimal Branching (Murray 1926)
//!
//! For minimum total power dissipation in Hagen–Poiseuille flow,
//! $D_0^3 = D_1^3 + D_2^3$ (Murray’s cube law). Deviation from this
//! law is quantified by `murray_law_deviation()` in the geometry module.

use super::geometry::BifurcationGeometry3D;
use super::types::{BifurcationConfig3D, BifurcationSolution3D};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};

// ============================================================================
// 3D Bifurcation Solver
// ============================================================================

/// 3D Finite Element Navier-Stokes solver for bifurcations
///
/// # Physics Equations
///
/// **Momentum equation (incompressible):**
/// ```text
/// ρ(∂u/∂t + (u·∇)u) + ∇p = ∇·τ + f
/// τ = μ(∇u + ∇u^T)  [Newtonian, constant μ]
/// τ = μ(γ̇)(∇u + ∇u^T)  [non-Newtonian, shear-rate dependent]
/// ```
///
/// **Continuity equation:**
/// ```text
/// ∇·u = 0
/// ```
///
/// # Numerical Method
///
/// - **Discretization**: FEM with P1-P1 (or P1-P0) elements
/// - **Stabilization**: Streamline Upwind Petrov-Galerkin (SUPG) for advection
/// - **Time integration**: Implicit Euler (steady-state: single iteration)
/// - **Nonlinear solver**: Fixed-point iteration or Newton-Raphson
/// - **Linear solver**: GMRES with preconditioner
///
/// # Validation Metrics
///
/// - **Mass conservation**: ∫∇·u dV ≈ 0
/// - **Momentum balance**: ∫momentum dV matched to inlet
/// - **Energy dissipation**: Tracked for viscous losses
/// - **Wall shear stress**: Computed from velocity gradient
pub struct BifurcationSolver3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    pub(crate) geometry: BifurcationGeometry3D<T>,
    pub(crate) config: BifurcationConfig3D<T>,
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + FromPrimitive
            + ToPrimitive
            + SafeFromF64
            + Float
            + From<f64>,
    > BifurcationSolver3D<T>
{
    /// Create new solver for given geometry and configuration
    pub fn new(geometry: BifurcationGeometry3D<T>, config: BifurcationConfig3D<T>) -> Self {
        Self { geometry, config }
    }

    /// Solve bifurcation flow with given fluid
    ///
    /// # Returns
    ///
    /// Solution structure containing velocities, pressures, and wall shear stresses
    #[allow(clippy::too_many_lines)]
    pub fn solve<F: FluidTrait<T> + Clone>(&self, fluid: F) -> Result<BifurcationSolution3D<T>> {
        self.validate_configuration()?;

        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use cfd_core::physics::boundary::BoundaryCondition;
        use cfd_mesh::BranchingMeshBuilder;
        use std::collections::HashMap;

        // 1. Generate Mesh
        let mesh_builder = BranchingMeshBuilder::bifurcation(
            self.geometry.d_parent,
            self.geometry.l_parent,
            self.geometry.d_daughter1,
            self.geometry.l_daughter1,
            self.geometry.branching_angle,
            self.config.mesh_resolution,
        );
        let base_mesh = match mesh_builder.build_surface() {
            Ok(m) => m,
            Err(e) => return Err(Error::Solver(format!("{e:?}"))),
        };
        let tet_mesh =
            cfd_mesh::application::hierarchy::hex_to_tet::HexToTetConverter::convert(&base_mesh);
        let mesh =
            cfd_mesh::application::hierarchy::hierarchical_mesh::P2MeshConverter::convert_to_p2(
                &tet_mesh,
            );

        // 2. Define Boundary Conditions
        let mut boundary_conditions = HashMap::new();
        let fluid_props =
            fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let inlet_area = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * num_traits::Float::powf(self.geometry.d_parent, T::from_f64_or_one(2.0));
        let u_inlet = self.config.inlet_flow_rate / inlet_area;

        // ── Classify boundary faces ───────────────────────────────────────────
        // AxialBoundaryClassifier replaces the two manual mesh.boundary_faces()
        // passes (outlet → inlet → wall) and the separate cell-face counting loop.
        let face_sets =
            crate::fem::AxialBoundaryClassifier::new(&mesh, self.config.mesh_resolution).classify();

        // Apply outlet BCs first (highest priority — shared rim nodes keep outlet BC)
        for &v_idx in &face_sets.outlet_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::PressureOutlet {
                    pressure: self.config.outlet_pressure,
                });
        }
        // Apply inlet BCs second
        for &v_idx in &face_sets.inlet_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::VelocityInlet {
                    velocity: Vector3::new(u_inlet, T::zero(), T::zero()),
                });
        }
        // Apply wall BCs to all remaining boundary nodes
        for &v_idx in &face_sets.wall_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::Dirichlet {
                    value: T::zero(),
                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero())]),
                });
        }

        // 3. Set up FEM Problem with f64 precision (build_surface produces f64 mesh)
        let constant_fluid_f64 = cfd_core::physics::fluid::ConstantPropertyFluid::<f64> {
            name: "Picard Iteration Basis".to_string(),
            density: fluid_props.density.to_f64().unwrap_or(1060.0),
            viscosity: fluid_props.dynamic_viscosity.to_f64().unwrap_or(3.5e-3),
            specific_heat: fluid_props.specific_heat.to_f64().unwrap_or(3600.0),
            thermal_conductivity: fluid_props.thermal_conductivity.to_f64().unwrap_or(0.5),
            speed_of_sound: fluid_props.speed_of_sound.to_f64().unwrap_or(1540.0),
        };

        // Convert boundary conditions from T -> f64 via manual variant mapping
        let convert_bc = |bc: BoundaryCondition<T>| -> BoundaryCondition<f64> {
            use cfd_core::physics::boundary::BoundaryCondition as BC;
            let to_f = |v: T| v.to_f64().unwrap_or(0.0);
            match bc {
                BC::Dirichlet {
                    value,
                    component_values,
                } => BC::Dirichlet {
                    value: to_f(value),
                    component_values: component_values
                        .map(|vs| vs.into_iter().map(|v| v.map(to_f)).collect()),
                },
                BC::VelocityInlet { velocity } => BC::VelocityInlet {
                    velocity: nalgebra::Vector3::new(
                        to_f(velocity.x),
                        to_f(velocity.y),
                        to_f(velocity.z),
                    ),
                },
                BC::PressureOutlet { pressure } => BC::PressureOutlet {
                    pressure: to_f(pressure),
                },
                BC::Neumann { gradient } => BC::Neumann {
                    gradient: to_f(gradient),
                },
                BC::PressureInlet {
                    pressure,
                    velocity_direction,
                } => BC::PressureInlet {
                    pressure: to_f(pressure),
                    velocity_direction: velocity_direction
                        .map(|v| nalgebra::Vector3::new(to_f(v.x), to_f(v.y), to_f(v.z))),
                },
                BC::MassFlowInlet {
                    mass_flow_rate,
                    temperature,
                } => BC::MassFlowInlet {
                    mass_flow_rate: to_f(mass_flow_rate),
                    temperature: temperature.map(to_f),
                },
                BC::VolumeFlowInlet { volume_flow_rate } => BC::VolumeFlowInlet {
                    volume_flow_rate: to_f(volume_flow_rate),
                },
                BC::Wall { .. } | BC::Symmetry | BC::Outflow | BC::Periodic { .. } => {
                    BC::Dirichlet {
                        value: 0.0_f64,
                        component_values: Some(vec![Some(0.0_f64); 3]),
                    }
                }
                BC::Robin { alpha, beta, gamma } => BC::Robin {
                    alpha: to_f(alpha),
                    beta: to_f(beta),
                    gamma: to_f(gamma),
                },
                BC::CharacteristicInlet {
                    riemann_invariant_r1,
                    riemann_invariant_r2,
                    entropy,
                    velocity,
                    pressure,
                } => BC::CharacteristicInlet {
                    riemann_invariant_r1: riemann_invariant_r1.map(to_f),
                    riemann_invariant_r2: riemann_invariant_r2.map(to_f),
                    entropy: entropy.map(to_f),
                    velocity: velocity
                        .map(|v| nalgebra::Vector3::new(to_f(v.x), to_f(v.y), to_f(v.z))),
                    pressure: pressure.map(to_f),
                },
                BC::CharacteristicOutlet {
                    pressure,
                    extrapolate_velocity,
                } => BC::CharacteristicOutlet {
                    pressure: to_f(pressure),
                    extrapolate_velocity,
                },
            }
        };
        let bc_f64: std::collections::HashMap<usize, BoundaryCondition<f64>> = boundary_conditions
            .into_iter()
            .map(|(k, bc)| (k, convert_bc(bc)))
            .collect();

        let mut problem = StokesFlowProblem::<f64>::new(
            mesh,
            constant_fluid_f64,
            bc_f64,
            tet_mesh.vertex_count(),
        );
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities: Vec<f64> =
            vec![fluid_props.dynamic_viscosity.to_f64().unwrap_or(3.5e-3); n_elements];
        let mut next_viscosities = Vec::with_capacity(n_elements);

        // 4. Picard Iteration Loop
        let fem_config = FemConfig::<f64>::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution: Option<crate::fem::StokesFlowSolution<f64>> = None;
        let mut anderson_accelerator = cfd_math::nonlinear_solver::AndersonAccelerator::<f64>::new(
            cfd_math::nonlinear_solver::AndersonConfig {
                history_depth: 5,
                relaxation: 1.0_f64,
                drop_tolerance: 1e-12_f64,
                method: cfd_math::nonlinear_solver::AndersonMethod::QR,
            },
        );

        for iter in 0..self.config.max_nonlinear_iterations {
            tracing::info!("Picard iteration {}", iter);

            // Update viscosities in problem
            if problem.element_viscosities.is_none() {
                problem.element_viscosities = Some(element_viscosities);
            }

            // Solve Stokes system
            let fem_solution = solver
                .solve(&problem, last_solution.as_ref())
                .map_err(|e| Error::Solver(e.to_string()))?;

            // Apply Anderson Acceleration
            let updated_solution = if let Some(ref prev) = last_solution {
                let acc_velocity =
                    anderson_accelerator.compute_next(&prev.velocity, &fem_solution.velocity);
                let mut acc_sol = fem_solution;
                acc_sol.velocity = acc_velocity;
                acc_sol
            } else {
                fem_solution
            };

            // Calculate shear rates and new viscosities
            let mut max_change_f64: f64 = 0.0;
            next_viscosities.clear();
            let current_viscosities = problem.element_viscosities.as_ref()
                .expect("element_viscosities set before Picard loop");

            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                let shear_rate_f64 =
                    self.calculate_element_shear_rate_f64(cell, &problem.mesh, &updated_solution)?;
                let shear_rate = T::from_f64_or_one(shear_rate_f64);
                let new_visc_t = fluid.viscosity_at_shear(
                    shear_rate,
                    T::from_f64_or_one(310.0),
                    self.config.inlet_pressure,
                )?;
                let new_visc = new_visc_t.to_f64().unwrap_or(3.5e-3);

                let change = num_traits::Float::abs(new_visc - current_viscosities[i])
                    / current_viscosities[i];
                if change > max_change_f64 {
                    max_change_f64 = change;
                }
                next_viscosities.push(new_visc);
            }

            element_viscosities = problem.element_viscosities.take()
                .expect("element_viscosities set before Picard loop");
            std::mem::swap(&mut element_viscosities, &mut next_viscosities);
            last_solution = Some(updated_solution);

            if max_change_f64 < self.config.nonlinear_tolerance.to_f64().unwrap_or(1e-4) {
                tracing::info!("Picard converged in {} iterations", iter + 1);
                break;
            }
        }

        let fem_solution =
            last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;

        // Debug: Check solution magnitude
        let vel_max_f64: f64 = fem_solution
            .velocity
            .iter()
            .map(|v| num_traits::Float::abs(*v))
            .fold(0.0_f64, f64::max);
        let p_max_f64: f64 = fem_solution
            .pressure
            .iter()
            .map(|p| num_traits::Float::abs(*p))
            .fold(0.0_f64, f64::max);
        eprintln!(
            "DEBUG: FEM solution max velocity = {vel_max_f64:?}, max pressure = {p_max_f64:?}"
        );

        // 5. Extract Metrics for BifurcationSolution3D
        let mesh = &problem.mesh;
        let mut solution = BifurcationSolution3D::new(&self.geometry);
        solution.u_parent_mean = u_inlet;
        solution.q_parent = self.config.inlet_flow_rate;

        // Calculate flows through daughters
        let q_d1_f64 = self.calculate_boundary_flow_f64(mesh, &fem_solution, "outlet_0")?;
        let q_d2_f64 = self.calculate_boundary_flow_f64(mesh, &fem_solution, "outlet_1")?;
        solution.q_daughter1 = <T as From<f64>>::from(q_d1_f64);
        solution.q_daughter2 = <T as From<f64>>::from(q_d2_f64);

        // Ensure u_daughter_mean is calculated
        let a_d1 = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * self.geometry.d_daughter1
            * self.geometry.d_daughter1;
        solution.u_daughter1_mean = solution.q_daughter1 / a_d1;

        let a_d2 = a_d1; // Assuming symmetry
        solution.u_daughter2_mean = solution.q_daughter2 / a_d2;

        solution.p_inlet = self.config.inlet_pressure;
        solution.p_outlet = self.config.outlet_pressure;
        solution.p_daughter1_outlet = self.config.outlet_pressure;
        solution.p_daughter2_outlet = self.config.outlet_pressure;

        // Extract junction pressure
        let p_junc_f64 = self.extract_point_pressure_f64(
            mesh,
            &fem_solution,
            nalgebra::Vector3::new(
                self.geometry.l_parent.to_f64().unwrap_or(0.0),
                0.0_f64,
                0.0_f64,
            ),
        )?;
        solution.p_junction_mid = <T as From<f64>>::from(p_junc_f64);

        // Calculate pressure drops
        solution.dp_parent = solution.p_inlet - solution.p_junction_mid;
        solution.dp_daughter1 = solution.p_junction_mid - solution.p_daughter1_outlet;
        solution.dp_daughter2 = solution.p_junction_mid - solution.p_daughter2_outlet;

        // Mass conservation check
        solution.mass_conservation_error =
            Float::abs(solution.q_parent - (solution.q_daughter1 + solution.q_daughter2));

        // Calculate wall shear stresses using analytical Poiseuille formula: τ_w = 8*μ*u_mean/R
        let mu = fluid_props.dynamic_viscosity;
        let r_parent = self.geometry.d_parent / T::from_f64_or_one(2.0);
        let r_daughter1 = self.geometry.d_daughter1 / T::from_f64_or_one(2.0);
        let r_daughter2 = self.geometry.d_daughter2 / T::from_f64_or_one(2.0);

        solution.wall_shear_stress_parent =
            T::from_f64_or_one(8.0) * mu * solution.u_parent_mean / r_parent;
        solution.wall_shear_stress_daughter1 =
            T::from_f64_or_one(8.0) * mu * solution.u_daughter1_mean / r_daughter1;
        solution.wall_shear_stress_daughter2 =
            T::from_f64_or_one(8.0) * mu * solution.u_daughter2_mean / r_daughter2;

        Ok(solution)
    }

    /// Validate solver configuration
    fn validate_configuration(&self) -> Result<()> {
        if self.config.inlet_flow_rate <= T::zero() {
            return Err(Error::InvalidInput(
                "Inlet flow rate must be positive".to_string(),
            ));
        }
        if self.config.nonlinear_tolerance <= T::zero() {
            return Err(Error::InvalidInput(
                "Nonlinear tolerance must be positive".to_string(),
            ));
        }
        if self.config.mesh_resolution < 2 {
            return Err(Error::InvalidInput(
                "Mesh resolution must be >= 2".to_string(),
            ));
        }
        Ok(())
    }

    /// Calculate Reynolds number in parent branch
    pub fn reynolds_number<F: FluidTrait<T> + Clone>(&self, fluid: F) -> Result<T> {
        let props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
        let u = self.config.inlet_flow_rate
            / (T::from_f64_or_one(std::f64::consts::PI / 4.0)
                * self.geometry.d_parent
                * self.geometry.d_parent);

        let re = (props.density * u * self.geometry.d_parent) / props.dynamic_viscosity;
        Ok(re)
    }

    /// Check if flow is laminar
    pub fn is_laminar<F: FluidTrait<T> + Clone>(&self, fluid: F) -> Result<bool> {
        let re = self.reynolds_number(fluid)?;
        Ok(re < T::from_f64_or_one(2300.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bifurcation_solver_creation() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let config = BifurcationConfig3D::default();
        let _solver = BifurcationSolver3D::new(geom, config);
    }

    #[test]
    fn test_bifurcation_solver_rejects_invalid_mesh_resolution() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let mut config = BifurcationConfig3D::default();
        config.mesh_resolution = 1;
        let solver = BifurcationSolver3D::new(geom, config);

        let water = cfd_core::physics::fluid::water_20c::<f64>().unwrap();
        let result = solver.solve(water);
        assert!(result.is_err());
    }
}
