//! 3D solver for flow in trifurcations
//!
//! Provides simulation for three-way branching vessels with non-Newtonian blood.
//!
//! # Theorem — Mass Conservation at Trifurcation
//!
//! For incompressible steady-state flow:
//!
//! ```text
//! Q_parent = Q_1 + Q_2 + Q_3
//! ```
//!
//! # Theorem — Generalised Murray’s Law ($N = 3$)
//!
//! $D_0^3 = D_1^3 + D_2^3 + D_3^3$ minimises total viscous dissipation
//! plus metabolic blood-volume maintenance cost.

use crate::trifurcation::geometry::TrifurcationGeometry3D;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::{Fluid as FluidTrait, NonNewtonianFluid};
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// Configuration for 3D trifurcation solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationConfig3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet volume flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Inlet pressure [Pa]
    pub inlet_pressure: T,
    /// Outlet pressures for each of the three daughter branches [Pa]
    pub outlet_pressures: [T; 3],
    /// Maximum Picard/Newton nonlinear iterations
    pub max_nonlinear_iterations: usize,
    /// Convergence tolerance for nonlinear residual
    pub nonlinear_tolerance: T,
    /// Maximum Krylov (GMRES/CG) linear iterations per solve
    pub max_linear_iterations: usize,
    /// Convergence tolerance for the linear solver
    pub linear_tolerance: T,
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + FromPrimitive
            + ToPrimitive
            + SafeFromF64,
    > Default for TrifurcationConfig3D<T>
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: T::from_f64_or_one(1e-8),
            inlet_pressure: T::from_f64_or_one(100.0),
            outlet_pressures: [T::zero(), T::zero(), T::zero()],
            max_nonlinear_iterations: 20,
            nonlinear_tolerance: T::from_f64_or_one(1e-4),
            max_linear_iterations: 1000,
            linear_tolerance: T::from_f64_or_one(1e-6),
        }
    }
}

/// Solution result for 3D trifurcation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationSolution3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Volume flow rates [parent, daughter1, daughter2, daughter3] [m³/s]
    pub flow_rates: [T; 4],
    /// Cross-section-averaged velocities [parent, d1, d2, d3] [m/s]
    pub mean_velocities: [T; 4],
    /// Wall shear stresses [parent, d1, d2, d3] [Pa]
    pub wall_shear_stresses: [T; 4],
    /// Pressure drops [parent, d1, d2, d3] [Pa]
    pub pressure_drops: [T; 4],
    /// Relative mass conservation error: |Q_in - ΣQ_out| / Q_in
    pub mass_conservation_error: T,
}

/// 3D Trifurcation flow solver
pub struct TrifurcationSolver3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Trifurcation channel geometry
    pub geometry: TrifurcationGeometry3D<T>,
    /// Solver configuration
    pub config: TrifurcationConfig3D<T>,
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
    > TrifurcationSolver3D<T>
{
    /// Create a new trifurcation solver with the given geometry and config
    pub fn new(geometry: TrifurcationGeometry3D<T>, config: TrifurcationConfig3D<T>) -> Self {
        Self { geometry, config }
    }

    /// Solve trifurcation flow
    #[allow(clippy::too_many_lines)]
    pub fn solve<F: FluidTrait<T> + NonNewtonianFluid<T>>(
        &self,
        fluid: &F,
    ) -> Result<TrifurcationSolution3D<T>> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use cfd_core::physics::boundary::BoundaryCondition;
        use cfd_mesh::BranchingMeshBuilder;
        use std::collections::HashMap;

        // 1. Generate Mesh
        let mesh_builder = BranchingMeshBuilder::trifurcation(
            self.geometry.d_parent,
            self.geometry.l_parent,
            self.geometry.d_daughters[0],
            self.geometry.l_daughters[0],
            self.geometry.branching_angles[0],
            8, // resolution factor
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

        let stats_vertex_count = mesh.vertex_count();
        let stats_cell_count = mesh.cell_count();
        let stats_boundary_face_count = mesh.boundary_faces().len();
        tracing::debug!(
            nodes = stats_vertex_count,
            cells = stats_cell_count,
            boundary_faces = stats_boundary_face_count,
            "Trifurcation mesh stats"
        );
        let mut label_counts = std::collections::HashMap::new();
        for f_idx in 0..mesh.face_count() {
            if let Some(label) = mesh.boundary_label(FaceId::from_usize(f_idx)) {
                *label_counts.entry(label.to_string()).or_insert(0) += 1;
            }
        }
        tracing::debug!(?label_counts, "Trifurcation boundary label counts");

        // 2. Define Boundary Conditions
        let mut boundary_conditions = HashMap::new();
        let fluid_props =
            fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let inlet_area = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * num_traits::Float::powf(self.geometry.d_parent, T::from_f64_or_one(2.0));
        let q_inlet_target = self.config.inlet_flow_rate;
        let u_inlet = q_inlet_target / inlet_area;

        // ── Classify boundary faces ───────────────────────────────────────────
        // AxialBoundaryClassifier replaces the three marked_faces iteration passes
        // (inlet / outlet_* / wall).  The per-label outlet map preserves per-outlet
        // pressure assignment needed for the trifurcation case.
        let face_sets = crate::fem::AxialBoundaryClassifier::new(&mesh, 8).classify();
        tracing::debug!(count = mesh.boundary_faces().len(), "Trifurcation boundary face count");

        // Apply inlet BCs (highest priority)
        for &v_idx in &face_sets.inlet_nodes {
            boundary_conditions.insert(
                v_idx,
                BoundaryCondition::VelocityInlet {
                    velocity: Vector3::new(u_inlet.to_f64().unwrap_or(0.0), 0.0_f64, 0.0_f64),
                },
            );
        }

        // Apply per-outlet pressure BCs using the per-label map (outlet_0, outlet_1, outlet_2)
        for (label, nodes) in &face_sets.outlet_nodes_by_label {
            // Parse the outlet index from the label suffix (e.g. "outlet_2" → 2).
            let idx = label
                .strip_prefix("outlet_")
                .and_then(|s| s.parse::<usize>().ok())
                .expect("Outlet label should have a valid integer suffix")
                .min(self.config.outlet_pressures.len().saturating_sub(1));
            let pressure = self.config.outlet_pressures[idx].to_f64().unwrap_or(0.0);
            for &v_idx in nodes {
                boundary_conditions
                    .entry(v_idx)
                    .or_insert(BoundaryCondition::PressureOutlet { pressure });
            }
        }

        // Apply wall (no-slip) BCs to remaining boundary nodes
        for &v_idx in &face_sets.wall_nodes {
            boundary_conditions
                .entry(v_idx)
                .or_insert(BoundaryCondition::Dirichlet {
                    value: 0.0_f64,
                    component_values: Some(vec![Some(0.0_f64), Some(0.0_f64), Some(0.0_f64), None]),
                });
        }

        tracing::debug!(count = boundary_conditions.len(), "Trifurcation boundary nodes constrained");

        // 3. Set up FEM Problem
        let constant_fluid = cfd_core::physics::fluid::ConstantPropertyFluid::<f64> {
            name: "Picard Iteration Basis".to_string(),
            density: fluid_props.density.to_f64().unwrap_or(0.0),
            viscosity: fluid_props.dynamic_viscosity.to_f64().unwrap_or(0.0),
            specific_heat: fluid_props.specific_heat.to_f64().unwrap_or(0.0),
            thermal_conductivity: fluid_props.thermal_conductivity.to_f64().unwrap_or(0.0),
            speed_of_sound: fluid_props.speed_of_sound.to_f64().unwrap_or(0.0),
        };

        let mut problem = StokesFlowProblem::<f64>::new(
            mesh,
            constant_fluid,
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

            let fem_solution = match fem_result {
                Ok(sol) => sol,
                Err(e) => {
                    tracing::warn!(error = %e, iter, "Trifurcation Picard linear solve failed; using last converged solution");
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
                    self.calculate_element_shear_rate_f64(cell, &problem.mesh, &updated_solution)?;
                let shear_rate = <T as From<f64>>::from(shear_rate_f64);
                let new_visc_t = fluid.apparent_viscosity(shear_rate);
                let new_visc = new_visc_t.to_f64().unwrap_or(0.0);
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
            tracing::debug!(iter, visc_change = max_change_f64, "Trifurcation Picard iteration");
            if max_change_f64 < self.config.nonlinear_tolerance.to_f64().unwrap_or(1e-4) {
                break;
            }
        }

        let fem_solution =
            last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        let mesh = &problem.mesh;

        // 5. Build Result
        let q_parent_fem = <T as From<f64>>::from(self.calculate_boundary_flow_f64(
            mesh,
            &fem_solution,
            "inlet",
        )?);
        let q_d1 = <T as From<f64>>::from(self.calculate_boundary_flow_f64(
            mesh,
            &fem_solution,
            "outlet_0",
        )?);
        let q_d2 = <T as From<f64>>::from(self.calculate_boundary_flow_f64(
            mesh,
            &fem_solution,
            "outlet_1",
        )?);
        let q_d3 = <T as From<f64>>::from(self.calculate_boundary_flow_f64(
            mesh,
            &fem_solution,
            "outlet_2",
        )?);

        let a_d1 = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * self.geometry.d_daughters[0]
            * self.geometry.d_daughters[0];
        let a_d2 = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * self.geometry.d_daughters[1]
            * self.geometry.d_daughters[1];
        let a_d3 = T::from_f64_or_one(std::f64::consts::PI / 4.0)
            * self.geometry.d_daughters[2]
            * self.geometry.d_daughters[2];

        let u_d1 = q_d1 / a_d1;
        let u_d2 = q_d2 / a_d2;
        let u_d3 = q_d3 / a_d3;

        let compute_poiseuille_wss = |q: T, d: T, mu: T| -> T {
            T::from_f64_or_one(32.0) * mu * q
                / (T::from_f64_or_one(std::f64::consts::PI) * d * d * d)
        };

        let mu_eff = fluid_props.dynamic_viscosity;
        let wss_parent = compute_poiseuille_wss(q_parent_fem, self.geometry.d_parent, mu_eff);
        let wss_d1 = compute_poiseuille_wss(q_d1, self.geometry.d_daughters[0], mu_eff);
        let wss_d2 = compute_poiseuille_wss(q_d2, self.geometry.d_daughters[1], mu_eff);
        let wss_d3 = compute_poiseuille_wss(q_d3, self.geometry.d_daughters[2], mu_eff);

        let wss = [wss_parent, wss_d1, wss_d2, wss_d3];

        let p_inlet_avg = {
            let mut sum = T::zero();
            let mut cnt: usize = 0;
            for f_idx in 0..mesh.face_count() {
                let fid = FaceId::from_usize(f_idx);
                if mesh.boundary_label(fid) == Some("inlet") {
                    let face = mesh.faces.get(fid);
                    for &vi in &face.vertices {
                        let vi_u = vi.as_usize();
                        if vi_u < fem_solution.n_corner_nodes {
                            sum += T::from_f64_or_one(fem_solution.get_pressure(vi_u));
                            cnt += 1;
                        }
                    }
                }
            }
            if cnt > 0 {
                sum / T::from_usize(cnt)
                    .expect("inlet pressure averaging count is always a representable usize")
            } else {
                self.config.inlet_pressure
            }
        };

        let extract_outlet_pressure = |label: &str| -> T {
            let mut sum = T::zero();
            let mut cnt: usize = 0;
            for f_idx in 0..mesh.face_count() {
                let fid = FaceId::from_usize(f_idx);
                if mesh.boundary_label(fid) == Some(label) {
                    let face = mesh.faces.get(fid);
                    for &vi in &face.vertices {
                        let vi_u = vi.as_usize();
                        if vi_u < fem_solution.n_corner_nodes {
                            sum += T::from_f64_or_one(fem_solution.get_pressure(vi_u));
                            cnt += 1;
                        }
                    }
                }
            }
            if cnt > 0 {
                sum / T::from_usize(cnt)
                    .expect("outlet pressure averaging count is always a representable usize")
            } else {
                T::zero()
            }
        };

        let p_out0 = extract_outlet_pressure("outlet_0");
        let p_out1 = extract_outlet_pressure("outlet_1");
        let p_out2 = extract_outlet_pressure("outlet_2");

        let p_out_mean = (p_out0 + p_out1 + p_out2) / T::from_f64_or_one(3.0);
        let dp = [
            Float::abs(p_inlet_avg - p_out_mean),
            Float::abs(p_inlet_avg - p_out0),
            Float::abs(p_inlet_avg - p_out1),
            Float::abs(p_inlet_avg - p_out2),
        ];

        let solution = TrifurcationSolution3D {
            flow_rates: [q_parent_fem, q_d1, q_d2, q_d3],
            mean_velocities: [u_inlet, u_d1, u_d2, u_d3],
            wall_shear_stresses: wss,
            pressure_drops: dp,
            mass_conservation_error: Float::abs(q_parent_fem - (q_d1 + q_d2 + q_d3)),
        };

        Ok(solution)
    }

    /// Calculate flow rate through a boundary label using u·n integration
    fn calculate_boundary_flow_f64(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        label: &str,
    ) -> Result<f64> {
        let reference_normal = self.boundary_reference_normal(label);
        let mut total_q = 0.0_f64;
        for f_idx in 0..mesh.face_count() {
            let fid = FaceId::from_usize(f_idx);
            if mesh.boundary_label(fid) == Some(label) {
                let face = mesh.faces.get(fid);
                if face.vertices.len() >= 3 {
                    let v0 = mesh.vertices.get(face.vertices[0]).position.coords;
                    let v1 = mesh.vertices.get(face.vertices[1]).position.coords;
                    let v2 = mesh.vertices.get(face.vertices[2]).position.coords;
                    let n_vec = (v1 - v0).cross(&(v2 - v0));
                    let area = n_vec.norm() * 0.5_f64;
                    if area <= 0.0_f64 {
                        continue;
                    }
                    let mut face_normal = n_vec.normalize();
                    if let Some(ref_n) = reference_normal {
                        if face_normal.dot(&ref_n) < 0.0_f64 {
                            face_normal = -face_normal;
                        }
                    }
                    let mut u_avg = nalgebra::Vector3::zeros();
                    for &v_idx in &face.vertices {
                        u_avg += solution.get_velocity(v_idx.as_usize());
                    }
                    u_avg /= face.vertices.len() as f64;
                    total_q += u_avg.dot(&face_normal) * area;
                }
            }
        }
        Ok(total_q.abs())
    }

    fn boundary_reference_normal(&self, label: &str) -> Option<nalgebra::Vector3<f64>> {
        match label {
            "inlet" => Some(nalgebra::Vector3::new(-1.0_f64, 0.0_f64, 0.0_f64)),
            "outlet_0" => {
                let theta = self.geometry.branching_angles[0].to_f64().unwrap_or(0.0);
                Some(nalgebra::Vector3::new(theta.cos(), theta.sin(), 0.0_f64))
            }
            "outlet_1" => {
                let theta = self.geometry.branching_angles[1].to_f64().unwrap_or(0.0);
                Some(nalgebra::Vector3::new(theta.cos(), theta.sin(), 0.0_f64))
            }
            "outlet_2" => {
                let theta = self.geometry.branching_angles[2].to_f64().unwrap_or(0.0);
                Some(nalgebra::Vector3::new(theta.cos(), theta.sin(), 0.0_f64))
            }
            _ => None,
        }
    }

    fn calculate_element_shear_rate_f64(
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
            use crate::fem::shape_functions::LagrangeTet10;
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
            let tet10 = LagrangeTet10::new(p1_grads);
            let l_centroid = [0.25_f64; 4];
            let p2_grads = tet10.gradients(&l_centroid);
            for i in 0..10 {
                let u = solution.get_velocity(idxs[i]);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += p2_grads[(col, i)] * u[row];
                    }
                }
            }
        } else {
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
        Ok((2.0_f64 * inner_prod).sqrt())
    }
}
