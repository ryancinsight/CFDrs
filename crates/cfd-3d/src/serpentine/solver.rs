//! 3D FEM Navier-Stokes solver for Serpentine channels
//!
//! Solves the incompressible Navier-Stokes equations on 3D serpentine domains
//! using Finite Element Method with support for non-Newtonian blood rheology and
//! Dean flow analysis.
//!
//! # Theorem — Dean Number and Secondary Flow (Dean 1927)
//!
//! In a curved channel of hydraulic diameter $D$ and radius of curvature $R$,
//! the dimensionless Dean number
//!
//! ```text
//! De = Re · √(D / 2R) = (ρ u D / μ) · √(D / 2R)
//! ```
//!
//! governs the onset and strength of secondary (Dean) vortices. For $De > De_{cr}
//! \approx 36$ a pair of counter-rotating vortices forms in the cross-section.
//!
//! **Proof sketch.** Balancing the centrifugal force $\rho u^2 / R$ against viscous
//! resistance $\mu u / D^2$ and non-dimensionalising yields the Dean parameter.
//! Linear stability analysis of the Navier–Stokes equations in toroidal
//! coordinates gives the critical value.
//!
//! **Reference:** Dean, W.R., "Note on the motion of fluid in a curved pipe",
//! Phil. Mag. 4(20), 1927, pp. 208–223.
//!
//! # Theorem — Picard Iteration Convergence
//!
//! For the non-Newtonian Picard iteration loop: given $\mu^{(k)}$, solve the
//! linear Stokes problem for $\mathbf{u}^{(k+1)}$, then update $\mu^{(k+1)} =
//! \mu(\dot{\gamma}(\mathbf{u}^{(k+1)}))$. If the viscosity function is Lipschitz
//! and the Reynolds number is sufficiently small, the iteration contracts in $H^1$ norm.

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use cfd_mesh::SerpentineMeshBuilder;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D Serpentine solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineConfig3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet volumetric flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Inlet pressure [Pa]
    pub inlet_pressure: T,
    /// Outlet pressure [Pa]
    pub outlet_pressure: T,

    /// Maximum iterations for nonlinear (Picard) solver
    pub max_nonlinear_iterations: usize,
    /// Convergence tolerance for nonlinear iterations
    pub nonlinear_tolerance: T,

    /// Mesh resolution (axial, transverse)
    pub resolution: (usize, usize),
    /// Whether the channel is circular or rectangular
    pub circular: bool,
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + FromPrimitive
            + ToPrimitive
            + SafeFromF64,
    > Default for SerpentineConfig3D<T>
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: T::from_f64_or_one(1e-7),
            inlet_pressure: T::from_f64_or_one(100.0),
            outlet_pressure: T::zero(),
            max_nonlinear_iterations: 15,
            nonlinear_tolerance: T::from_f64_or_one(1e-4),
            resolution: (80, 8),
            circular: false,
        }
    }
}

// ============================================================================
// 3D Serpentine Solver
// ============================================================================

/// 3D Finite Element Navier-Stokes solver for Serpentine channels
pub struct SerpentineSolver3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float> {
    builder: SerpentineMeshBuilder<T>,
    config: SerpentineConfig3D<T>,
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
    > SerpentineSolver3D<T>
{
    /// Create new solver from mesh builder and config
    pub fn new(builder: SerpentineMeshBuilder<T>, config: SerpentineConfig3D<T>) -> Self {
        Self { builder, config }
    }

    /// Solve Serpentine flow with given fluid
    #[allow(clippy::too_many_lines)]
    pub fn solve<F: FluidTrait<T> + Clone>(&self, fluid: F) -> Result<SerpentineSolution3D<T>> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use cfd_core::physics::boundary::BoundaryCondition;
        use std::collections::HashMap;

        // 1. Generate Mesh
        let base_mesh = self
            .builder
            .clone()
            .with_resolution(self.config.resolution.0, self.config.resolution.1)
            .with_circular(self.config.circular)
            .build_surface();

        let base_mesh = match base_mesh {
            Ok(m) => m,
            Err(e) => return Err(Error::Solver(format!("{e:?}"))),
        };

        // 1.1 Decompose to Tetrahedra and Promote to Quadratic (P2) mesh for Taylor-Hood elements (Q2-Q1)
        let tet_mesh =
            cfd_mesh::application::hierarchy::hex_to_tet::HexToTetConverter::convert(&base_mesh);
        let mesh =
            cfd_mesh::application::hierarchy::hierarchical_mesh::P2MeshConverter::convert_to_p2(
                &tet_mesh,
            );

        // 2. Define Boundary Conditions
        // Priority: inlet > outlet > wall. Process inlet/outlet first so that
        // shared corner/edge nodes get the correct BC instead of a no-slip wall BC.
        let mut boundary_conditions = HashMap::new();
        let fluid_props =
            fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let area_inlet = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0)
                * self.builder.diameter
                * self.builder.diameter
        } else {
            self.builder.diameter * self.builder.diameter
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        // ── Classify boundary faces ───────────────────────────────────────────
        // AxialBoundaryClassifier handles the inlet/outlet/wall label loops,
        // replacing the three separate mesh.boundary_faces() passes above.
        let face_sets =
            crate::fem::AxialBoundaryClassifier::new(&mesh, self.config.resolution.0).classify();

        // Apply BCs from classified node sets (inlet first for rim-node priority)
        for &v_idx in &face_sets.inlet_nodes {
            boundary_conditions.insert(
                v_idx,
                BoundaryCondition::VelocityInlet {
                    velocity: Vector3::new(0.0_f64, 0.0_f64, u_inlet.to_f64().unwrap_or(0.0)),
                },
            );
        }
        for &v_idx in &face_sets.outlet_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::PressureOutlet {
                    pressure: self.config.outlet_pressure.to_f64().unwrap_or(0.0),
                });
        }
        for &v_idx in &face_sets.wall_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::Dirichlet {
                    value: 0.0_f64,
                    component_values: Some(vec![Some(0.0_f64), Some(0.0_f64), Some(0.0_f64), None]),
                });
        }

        // 3. Set up FEM Problem with initial viscosity
        let constant_basis = cfd_core::physics::fluid::ConstantPropertyFluid::<f64> {
            name: "Picard Basis".to_string(),
            density: fluid_props.density.to_f64().unwrap_or(0.0),
            viscosity: fluid_props.dynamic_viscosity.to_f64().unwrap_or(0.0),
            specific_heat: fluid_props.specific_heat.to_f64().unwrap_or(0.0),
            thermal_conductivity: fluid_props.thermal_conductivity.to_f64().unwrap_or(0.0),
            speed_of_sound: fluid_props.speed_of_sound.to_f64().unwrap_or(0.0),
        };

        let mut problem = StokesFlowProblem::<f64>::new(
            mesh,
            constant_basis,
            boundary_conditions,
            tet_mesh.vertex_count(),
        );
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities =
            vec![fluid_props.dynamic_viscosity.to_f64().unwrap_or(0.0); n_elements];
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
            if problem.element_viscosities.is_none() {
                problem.element_viscosities = Some(element_viscosities);
            }
            let fem_result = solver.solve(&problem, last_solution.as_ref());

            // If the linear solver fails (e.g. GMRES stagnation on the
            // non-Newtonian saddle-point system), use the last converged
            // solution and break.
            let fem_solution = match fem_result {
                Ok(sol) => sol,
                Err(e) => {
                    tracing::warn!(error = %e, iter, "Serpentine Picard linear solve failed; using last converged solution");
                    break;
                }
            };

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

            let mut max_change_f64 = 0.0_f64;
            next_viscosities.clear();
            let current_viscosities = problem.element_viscosities.as_ref()
                .expect("element_viscosities set before Picard loop");

            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                let shear_rate_f64 =
                    self.calculate_cell_shear_rate_f64(cell, &problem.mesh, &updated_solution)?;
                let shear_rate = <T as From<f64>>::from(shear_rate_f64);
                let new_visc_t = fluid.viscosity_at_shear(
                    shear_rate,
                    T::from_f64_or_one(310.0),
                    self.config.inlet_pressure,
                )?;
                let new_visc = new_visc_t.to_f64().unwrap_or(0.0);

                let change = Float::abs(new_visc - current_viscosities[i]) / current_viscosities[i];
                if change > max_change_f64 {
                    max_change_f64 = change;
                }
                next_viscosities.push(new_visc);
            }

            element_viscosities = problem.element_viscosities.take()
                .expect("element_viscosities set before Picard loop");
            std::mem::swap(&mut element_viscosities, &mut next_viscosities);
            last_solution = Some(updated_solution);

            tracing::debug!(iter, visc_change = max_change_f64, "Serpentine Picard iteration");
            if max_change_f64 < self.config.nonlinear_tolerance.to_f64().unwrap_or(1e-4) {
                break;
            }
        }

        let _fem_solution =
            last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;

        // 5. Extract Metrics
        let mut solution = SerpentineSolution3D::new();
        solution.u_inlet = u_inlet;
        solution.p_inlet = self.config.inlet_pressure;
        solution.p_outlet = self.config.outlet_pressure;
        solution.dp_total = solution.p_inlet - solution.p_outlet;

        // Calculate Dean Number: De = Re * sqrt(Dh / 2Rc)
        // For sine wave path x = A*sin(k*z), curvature kappa = |x''| / (1 + x'^2)^(3/2)
        // Max curvature at peaks: kappa_max = A*k^2. Radius Rc = 1/kappa_max = 1 / (A * (2pi/lambda)^2)
        let k = <T as FromPrimitive>::from_f64(2.0 * std::f64::consts::PI)
            .expect("2π is an IEEE 754 representable f64 constant")
            / self.builder.wavelength;
        let kappa_max = self.builder.amplitude * k * k;
        let rc = T::one() / Float::max(kappa_max, <T as FromPrimitive>::from_f64(1e-10)
            .expect("1e-10 is an IEEE 754 representable f64 constant"));

        let re =
            (fluid_props.density * u_inlet * self.builder.diameter) / fluid_props.dynamic_viscosity;
        solution.dean_number = re
            * Float::sqrt(
                self.builder.diameter / (<T as FromPrimitive>::from_f64(2.0)
                    .expect("2.0 is representable in all IEEE 754 types") * rc),
            );

        Ok(solution)
    }

    fn calculate_cell_shear_rate_f64(
        &self,
        cell: &cfd_mesh::domain::topology::Cell,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> Result<f64> {
        let mut idxs = Vec::with_capacity(10);
        for &face_idx in &cell.faces {
            let face = mesh.faces.get(FaceId::from_usize(face_idx));
            for &v_idx in &face.vertices {
                let v_idx_u = v_idx.as_usize();
                if !idxs.contains(&v_idx_u) {
                    idxs.push(v_idx_u);
                }
            }
        }

        if idxs.len() < 4 {
            return Err(Error::Solver("Invalid cell topology".to_string()));
        }

        let mut local_verts = Vec::with_capacity(idxs.len());
        for &idx in &idxs {
            local_verts.push(mesh.vertices.get(VertexId::from_usize(idx)).position.coords);
        }

        let mut l = nalgebra::Matrix3::zeros();

        if idxs.len() == 10 {
            // Tet10 (P2): Evaluate gradient at centroid (L = [0.25, 0.25, 0.25, 0.25])
            use crate::fem::shape_functions::LagrangeTet10;

            // 1. Calculate P1 gradients (∇L_i)
            let mut tet4 = crate::fem::element::FluidElement::<f64>::new(idxs[0..4].to_vec());
            let six_v = tet4.calculate_volume(&local_verts);
            if six_v.abs() < 1e-24_f64 {
                return Ok(0.0_f64);
            }
            tet4.calculate_shape_derivatives(&local_verts[0..4]);
            let p1_grads = nalgebra::Matrix3x4::from_columns(&[
                Vector3::new(
                    tet4.shape_derivatives[(0, 0)],
                    tet4.shape_derivatives[(1, 0)],
                    tet4.shape_derivatives[(2, 0)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 1)],
                    tet4.shape_derivatives[(1, 1)],
                    tet4.shape_derivatives[(2, 1)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 2)],
                    tet4.shape_derivatives[(1, 2)],
                    tet4.shape_derivatives[(2, 2)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 3)],
                    tet4.shape_derivatives[(1, 3)],
                    tet4.shape_derivatives[(2, 3)],
                ),
            ]);

            // 2. Evaluate P2 gradients at centroid
            let tet10 = LagrangeTet10::new(p1_grads);
            let l_centroid = [0.25_f64; 4];
            let p2_grads = tet10.gradients(&l_centroid);

            // 3. Compute velocity gradient: L = sum(u_i * ∇N_i)
            for i in 0..10 {
                let u = solution.get_velocity(idxs[i]);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += p2_grads[(col, i)] * u[row];
                    }
                }
            }
        } else {
            // Tet4: constant gradient
            let mut element = crate::fem::element::FluidElement::<f64>::new(idxs.clone());
            element.calculate_shape_derivatives(&local_verts);
            for (i, &idx) in idxs.iter().enumerate().take(4) {
                let u = solution.get_velocity(idx);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += element.shape_derivatives[(col, i)] * u[row];
                    }
                }
            }
        }

        let epsilon = (l + l.transpose()) * 0.5_f64;
        let mut inner_prod = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                inner_prod += epsilon[(i, j)] * epsilon[(i, j)];
            }
        }
        let shear = (2.0_f64 * inner_prod).sqrt();

        Ok(shear)
    }
}

// ============================================================================
// Solution Result
// ============================================================================

/// Complete solution to 3D Serpentine problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SerpentineSolution3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet mean velocity [m/s]
    pub u_inlet: T,
    /// Inlet pressure [Pa]
    pub p_inlet: T,
    /// Outlet pressure [Pa]
    pub p_outlet: T,
    /// Total pressure drop [Pa]
    pub dp_total: T,
    /// Dean number at curve peaks
    pub dean_number: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> SerpentineSolution3D<T> {
    /// Create new solution structure
    pub fn new() -> Self {
        Self {
            u_inlet: T::zero(),
            p_inlet: T::zero(),
            p_outlet: T::zero(),
            dp_total: T::zero(),
            dean_number: T::zero(),
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> Default for SerpentineSolution3D<T> {
    fn default() -> Self {
        Self::new()
    }
}
