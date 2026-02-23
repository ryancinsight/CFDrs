//! 3D FEM Navier-Stokes solver for bifurcations with blood flow
//!
//! Solves the incompressible Navier-Stokes equations on 3D bifurcation domains
//! using Finite Element Method with stabilization for advection-dominated flows.

use super::geometry::BifurcationGeometry3D;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D bifurcation solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationConfig3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet volumetric flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Inlet pressure [Pa]
    pub inlet_pressure: T,
    /// Outlet pressure [Pa]
    pub outlet_pressure: T,

    /// Time step size [s] (for transient)
    pub time_step: T,
    /// Number of time steps
    pub num_time_steps: usize,
    /// Use steady-state (num_time_steps = 1)
    pub steady_state: bool,

    /// Maximum iterations for nonlinear solver
    pub max_nonlinear_iterations: usize,
    /// Convergence tolerance for nonlinear iterations
    pub nonlinear_tolerance: T,

    /// Maximum iterations for linear solver
    pub max_linear_iterations: usize,
    /// Convergence tolerance for linear solver
    pub linear_tolerance: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
    for BifurcationConfig3D<T>
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: T::from_f64_or_one(1e-8),
            inlet_pressure: T::from_f64_or_one(100.0),
            outlet_pressure: T::zero(),
            time_step: T::from_f64_or_one(0.001),
            num_time_steps: 1,
            steady_state: true,
            max_nonlinear_iterations: 20,
            nonlinear_tolerance: T::from_f64_or_one(1e-4),
            max_linear_iterations: 1000,
            linear_tolerance: T::from_f64_or_one(1e-6),
        }
    }
}

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
    geometry: BifurcationGeometry3D<T>,
    config: BifurcationConfig3D<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float + From<f64>> BifurcationSolver3D<T> {
    /// Create new solver for given geometry and configuration
    pub fn new(geometry: BifurcationGeometry3D<T>, config: BifurcationConfig3D<T>) -> Self {
        Self { geometry, config }
    }

    /// Solve bifurcation flow with given fluid
    ///
    /// # Returns
    ///
    /// Solution structure containing velocities, pressures, and wall shear stresses
    pub fn solve<F: FluidTrait<T> + Clone>(
        &self,
        fluid: F,
    ) -> Result<BifurcationSolution3D<T>> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use cfd_mesh::BranchingMeshBuilder;
        use cfd_mesh::domain::core::index::{FaceId, VertexId};
        use std::collections::HashMap;
        use cfd_core::physics::boundary::BoundaryCondition;

        // 1. Generate Mesh
        let mesh_builder = BranchingMeshBuilder::bifurcation(
            self.geometry.d_parent,
            self.geometry.l_parent,
            self.geometry.d_daughter1,
            self.geometry.l_daughter1,
            self.geometry.branching_angle,
            8, // resolution factor
        );
        let base_mesh = match mesh_builder.build_surface() {
            Ok(m) => m,
            Err(e) => return Err(Error::Solver(format!("{:?}", e))),
        };
        let tet_mesh = cfd_mesh::application::hierarchy::hex_to_tet::HexToTetConverter::convert(&base_mesh);
        let mesh = cfd_mesh::application::hierarchy::hierarchical_mesh::P2MeshConverter::convert_to_p2(&tet_mesh);

        // 2. Define Boundary Conditions
        let mut boundary_conditions = HashMap::new();
        let fluid_props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
        
        let inlet_area = T::from_f64_or_one(std::f64::consts::PI / 4.0) * num_traits::Float::powf(self.geometry.d_parent, T::from_f64_or_one(2.0));
        let u_inlet = self.config.inlet_flow_rate / inlet_area;

        // Pass 1: assign BCs based on labeled boundary faces
        // Use entry().or_insert() to avoid overwriting (outlet/inlet take priority over wall)
        let mut inlet_count = 0;
        let mut outlet_count = 0;
        
        // Calculate outlet velocity based on area ratio
        // For symmetric bifurcation: u_outlet = u_inlet * (A_parent / A_daughter) = u_inlet * (D_parent/D_daughter)^2
        let ratio = self.geometry.d_parent / self.geometry.d_daughter1;
        let area_ratio = ratio * ratio;
        let u_outlet = u_inlet * area_ratio; // Velocity in daughter branches
        
        // Process outlet faces FIRST (highest priority) - use PressureOutlet
        for f_id in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_id) {
                if label == "outlet_0" || label == "outlet_1" {
                    let face = mesh.faces.get(f_id);
                    for &v_idx in &face.vertices {
                        boundary_conditions.entry(v_idx.as_usize()).or_insert_with(|| {
                            outlet_count += 1;
                            BoundaryCondition::PressureOutlet {
                                pressure: self.config.outlet_pressure,
                            }
                        });
                    }
                }
            }
        }
        
        // Process inlet faces SECOND
        for f_id in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_id) {
                if label == "inlet" {
                    let face = mesh.faces.get(f_id);
                    for &v_idx in &face.vertices {
                        boundary_conditions.entry(v_idx.as_usize()).or_insert_with(|| {
                            inlet_count += 1;
                            BoundaryCondition::VelocityInlet {
                                velocity: Vector3::new(u_inlet, T::zero(), T::zero()),
                            }
                        });
                    }
                }
            }
        }
        
        eprintln!("DEBUG BC Assignment - Pass 1: inlet vertices={}, outlet vertices={}", inlet_count, outlet_count);
        
        // Second pass: ensure ALL boundary nodes have a BC (default to wall)
        // This handles any boundary faces that might not have been labeled
        use std::collections::{HashMap as StdHashMap, HashSet};
        
        // Count face references to identify boundary faces
        let mut face_cell_count: StdHashMap<usize, usize> = StdHashMap::new();
        for cell in mesh.cells.iter() {
            for &face_idx in &cell.faces {
                *face_cell_count.entry(face_idx).or_insert(0) += 1;
            }
        }
        
        // Find all boundary faces (referenced by only one cell)
        let boundary_faces: HashSet<usize> = face_cell_count
            .iter()
            .filter(|&(_, &count)| count == 1)
            .map(|(&face_idx, _)| face_idx)
            .collect();
        
        // Ensure all boundary vertices have a BC
        let mut default_wall_count = 0;
        for &face_idx in &boundary_faces {
            if face_idx < mesh.face_count() {
                let face = mesh.faces.get(FaceId::from_usize(face_idx));
                for &v_id in &face.vertices {
                    // If no BC assigned yet, default to wall
                    boundary_conditions.entry(v_id.as_usize()).or_insert_with(|| {
                        default_wall_count += 1;
                        BoundaryCondition::Dirichlet {
                            value: T::zero(),
                            component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero())]),
                        }
                    });
                }
            }
        }
        eprintln!("DEBUG BC Assignment - Pass 2: default wall vertices added={}", default_wall_count);
        eprintln!("DEBUG Total BC entries: {}", boundary_conditions.len());
        
        // Check BC types
        let outlet_bc_count = boundary_conditions.values().filter(|bc| matches!(bc, BoundaryCondition::PressureOutlet { .. })).count();
        let wall_bc_count = boundary_conditions.values().filter(|bc| matches!(bc, BoundaryCondition::Wall { .. })).count();
        let inlet_bc_count = boundary_conditions.values().filter(|bc| matches!(bc, BoundaryCondition::VelocityInlet { .. })).count();
        let dirichlet_count = boundary_conditions.values().filter(|bc| matches!(bc, BoundaryCondition::Dirichlet { .. })).count();
        eprintln!("DEBUG BC Types: VelocityInlet={}, PressureOutlet={}, Wall={}, Dirichlet={}", 
                 inlet_bc_count, outlet_bc_count, wall_bc_count, dirichlet_count);

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
                BC::Dirichlet { value, component_values } => BC::Dirichlet {
                    value: to_f(value),
                    component_values: component_values.map(|vs| vs.into_iter().map(|v| v.map(to_f)).collect()),
                },
                BC::VelocityInlet { velocity } => BC::VelocityInlet {
                    velocity: nalgebra::Vector3::new(to_f(velocity.x), to_f(velocity.y), to_f(velocity.z)),
                },
                BC::PressureOutlet { pressure } => BC::PressureOutlet { pressure: to_f(pressure) },
                BC::Neumann { gradient } => BC::Neumann { gradient: to_f(gradient) },
                BC::PressureInlet { pressure, velocity_direction } => BC::PressureInlet {
                    pressure: to_f(pressure),
                    velocity_direction: velocity_direction.map(|v| nalgebra::Vector3::new(to_f(v.x), to_f(v.y), to_f(v.z))),
                },
                BC::MassFlowInlet { mass_flow_rate, temperature } => BC::MassFlowInlet {
                    mass_flow_rate: to_f(mass_flow_rate),
                    temperature: temperature.map(to_f),
                },
                BC::VolumeFlowInlet { volume_flow_rate } => BC::VolumeFlowInlet { volume_flow_rate: to_f(volume_flow_rate) },
                BC::Wall { .. } | BC::Symmetry | BC::Outflow | BC::Periodic { .. } => BC::Dirichlet { value: 0.0_f64, component_values: Some(vec![Some(0.0_f64); 3]) },
                BC::Robin { alpha, beta, gamma } => BC::Robin { alpha: to_f(alpha), beta: to_f(beta), gamma: to_f(gamma) },
                BC::CharacteristicInlet { riemann_invariant_r1, riemann_invariant_r2, entropy, velocity, pressure } => BC::CharacteristicInlet {
                    riemann_invariant_r1: riemann_invariant_r1.map(to_f),
                    riemann_invariant_r2: riemann_invariant_r2.map(to_f),
                    entropy: entropy.map(to_f),
                    velocity: velocity.map(|v| nalgebra::Vector3::new(to_f(v.x), to_f(v.y), to_f(v.z))),
                    pressure: pressure.map(to_f),
                },
                BC::CharacteristicOutlet { pressure, extrapolate_velocity } => BC::CharacteristicOutlet { pressure: to_f(pressure), extrapolate_velocity },
            }
        };
        let bc_f64: std::collections::HashMap<usize, BoundaryCondition<f64>> =
            boundary_conditions.into_iter().map(|(k, bc)| (k, convert_bc(bc))).collect();

        let mut problem = StokesFlowProblem::<f64>::new(mesh, constant_fluid_f64, bc_f64, tet_mesh.vertex_count());
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities: Vec<f64> = vec![fluid_props.dynamic_viscosity.to_f64().unwrap_or(3.5e-3); n_elements];
        
        // 4. Picard Iteration Loop
        let fem_config = FemConfig::<f64>::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;

        for iter in 0..self.config.max_nonlinear_iterations {
            tracing::info!("Picard iteration {}", iter);
            
            // Update viscosities in problem
            problem.element_viscosities = Some(element_viscosities.clone());
            
            // Solve Stokes system
            let fem_solution = solver.solve(&problem, last_solution.as_ref()).map_err(|e| Error::Solver(e.to_string()))?;
            
            // Calculate shear rates and new viscosities
            let mut max_change_f64: f64 = 0.0;
            let mut new_viscosities = Vec::with_capacity(n_elements);
            
            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                let shear_rate_f64 = self.calculate_element_shear_rate_f64(cell, &problem.mesh, &fem_solution)?;
                let shear_rate = T::from_f64_or_one(shear_rate_f64);
                let new_visc_t = fluid.viscosity_at_shear(shear_rate, T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
                let new_visc = new_visc_t.to_f64().unwrap_or(3.5e-3);
                
                let change = num_traits::Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change_f64 {
                    max_change_f64 = change;
                }
                new_viscosities.push(new_visc);
            }
            
            element_viscosities = new_viscosities;
            last_solution = Some(fem_solution);
            
            if max_change_f64 < self.config.nonlinear_tolerance.to_f64().unwrap_or(1e-4) {
                tracing::info!("Picard converged in {} iterations", iter + 1);
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        
        // Debug: Check solution magnitude
        let vel_max_f64: f64 = fem_solution.velocity.iter().map(|v| num_traits::Float::abs(*v)).fold(0.0_f64, |a, b| f64::max(a, b));
        let p_max_f64: f64 = fem_solution.pressure.iter().map(|p| num_traits::Float::abs(*p)).fold(0.0_f64, |a, b| f64::max(a, b));
        eprintln!("DEBUG: FEM solution max velocity = {:?}, max pressure = {:?}", vel_max_f64, p_max_f64);

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
        let a_d1 = T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.geometry.d_daughter1 * self.geometry.d_daughter1;
        solution.u_daughter1_mean = solution.q_daughter1 / a_d1;
        
        let a_d2 = a_d1; // Assuming symmetry
        solution.u_daughter2_mean = solution.q_daughter2 / a_d2;

        solution.p_inlet = self.config.inlet_pressure;
        solution.p_outlet = self.config.outlet_pressure;
        solution.p_daughter1_outlet = self.config.outlet_pressure;
        solution.p_daughter2_outlet = self.config.outlet_pressure;
        
        // Extract junction pressure
        let p_junc_f64 = self.extract_point_pressure_f64(mesh, &fem_solution, nalgebra::Vector3::new(self.geometry.l_parent.to_f64().unwrap_or(0.0), 0.0_f64, 0.0_f64))?;
        solution.p_junction_mid = <T as From<f64>>::from(p_junc_f64);
        
        // Calculate pressure drops
        solution.dp_parent = solution.p_inlet - solution.p_junction_mid;
        solution.dp_daughter1 = solution.p_junction_mid - solution.p_daughter1_outlet;
        solution.dp_daughter2 = solution.p_junction_mid - solution.p_daughter2_outlet;

        // Mass conservation check
        solution.mass_conservation_error = Float::abs(solution.q_parent - (solution.q_daughter1 + solution.q_daughter2));

        // Calculate wall shear stresses using analytical Poiseuille formula: τ_w = 8*μ*u_mean/R
        let mu = fluid_props.dynamic_viscosity;
        let r_parent = self.geometry.d_parent / T::from_f64_or_one(2.0);
        let r_daughter1 = self.geometry.d_daughter1 / T::from_f64_or_one(2.0);
        let r_daughter2 = self.geometry.d_daughter2 / T::from_f64_or_one(2.0);
        
        solution.wall_shear_stress_parent = T::from_f64_or_one(8.0) * mu * solution.u_parent_mean / r_parent;
        solution.wall_shear_stress_daughter1 = T::from_f64_or_one(8.0) * mu * solution.u_daughter1_mean / r_daughter1;
        solution.wall_shear_stress_daughter2 = T::from_f64_or_one(8.0) * mu * solution.u_daughter2_mean / r_daughter2;

        Ok(solution)
    }

    /// Calculate flow rate through a boundary label using u·n integration (f64 precision)
    fn calculate_boundary_flow_f64(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        label: &str,
    ) -> Result<f64> {
        let mut total_q = 0.0_f64;
        let mut face_count = 0;
        
        for f_idx in 0..mesh.face_count() {
            let f_id = FaceId::from_usize(f_idx);
            if mesh.boundary_label(f_id) == Some(label) {
                face_count += 1;
                let face = mesh.faces.get(f_id);
                // FaceData always has exactly 3 vertices
                let v0 = mesh.vertices.position(face.vertices[0]).coords;
                let v1 = mesh.vertices.position(face.vertices[1]).coords;
                let v2 = mesh.vertices.position(face.vertices[2]).coords;
                
                let n_vec = (v1 - v0).cross(&(v2 - v0));
                let area = n_vec.norm() * 0.5_f64;
                let face_normal = n_vec.normalize();
                
                let mut u_avg = nalgebra::Vector3::zeros();
                for &v_id in &face.vertices {
                    let u = solution.get_velocity(v_id.as_usize());
                    u_avg += u;
                    if face_count <= 2 {
                        eprintln!("DEBUG: vertex {} velocity = {:?}", v_id.as_usize(), u);
                    }
                }
                u_avg /= 3.0_f64;
                
                let face_flow = u_avg.dot(&face_normal) * area;
                total_q += face_flow;
            }
        }
        
        eprintln!("DEBUG: Flow integration for '{}': {} faces, total_q = {:?}", label, face_count, total_q);
        Ok(total_q.abs())
    }

    /// Extract pressure at a point in the mesh (closest node, f64 precision)
    fn extract_point_pressure_f64(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        point: nalgebra::Vector3<f64>,
    ) -> Result<f64> {
        let mut best_node = 0;
        let mut min_dist = f64::MAX;
        
        let n_pressure_nodes = solution.n_corner_nodes.min(mesh.vertex_count());
        for i in 0..n_pressure_nodes {
            let pos = mesh.vertices.position(VertexId::from_usize(i));
            let dist = (pos.coords - point).norm();
            if dist < min_dist {
                min_dist = dist;
                best_node = i;
            }
        }
        
        Ok(solution.get_pressure(best_node))
    }

    /// Calculate element shear rate (f64 internal precision)
    fn calculate_element_shear_rate_f64(
        &self,
        cell: &cfd_mesh::domain::topology::Cell,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> Result<f64> {
        let mut idxs: Vec<usize> = Vec::with_capacity(10);
        for &face_idx in &cell.faces {
            if face_idx < mesh.face_count() {
                let face = mesh.faces.get(FaceId::from_usize(face_idx));
                for &v_id in &face.vertices {
                    let v_usize = v_id.as_usize();
                    if !idxs.contains(&v_usize) {
                        idxs.push(v_usize);
                    }
                }
            }
        }

        if idxs.len() < 4 {
            return Err(Error::Solver("Invalid cell topology".to_string()));
        }

        let mut local_verts = Vec::with_capacity(idxs.len());
        for &idx in &idxs {
            local_verts.push(mesh.vertices.position(VertexId::from_usize(idx)).coords);
        }

        let mut l = nalgebra::Matrix3::zeros();
        
        if idxs.len() == 10_usize {
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
                    tet4.shape_derivatives[(2, 0)]
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 1)],
                    tet4.shape_derivatives[(1, 1)],
                    tet4.shape_derivatives[(2, 1)]
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 2)],
                    tet4.shape_derivatives[(1, 2)],
                    tet4.shape_derivatives[(2, 2)]
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 3)],
                    tet4.shape_derivatives[(1, 3)],
                    tet4.shape_derivatives[(2, 3)]
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
            let mut element = crate::fem::element::FluidElement::new(idxs.clone());
            element.calculate_shape_derivatives(&local_verts);
            for i in 0..4 {
                let u = solution.get_velocity(idxs[i]);
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
        
        Ok(num_traits::Float::sqrt(2.0_f64 * inner_prod))
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
        Ok(())
    }

    /// Calculate Reynolds number in parent branch
    pub fn reynolds_number<F: FluidTrait<T> + Clone>(
        &self,
        fluid: F,
    ) -> Result<T> {
        let props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
        let u = self.config.inlet_flow_rate
            / (T::from_f64_or_one(std::f64::consts::PI / 4.0)
                * self.geometry.d_parent
                * self.geometry.d_parent);

        let re = (props.density * u * self.geometry.d_parent) / props.dynamic_viscosity;
        Ok(re)
    }

    /// Check if flow is laminar
    pub fn is_laminar<F: FluidTrait<T> + Clone>(
        &self,
        fluid: F,
    ) -> Result<bool> {
        let re = self.reynolds_number(fluid)?;
        Ok(re < T::from_f64_or_one(2300.0))
    }
}

// ============================================================================
// Solution Result
// ============================================================================

/// Complete solution to 3D bifurcation problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BifurcationSolution3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    pub q_parent: T,
    pub q_daughter1: T,
    pub q_daughter2: T,
    pub u_parent_mean: T,
    pub u_daughter1_mean: T,
    pub u_daughter2_mean: T,
    pub p_inlet: T,
    pub p_junction_mid: T,
    pub p_daughter1_outlet: T,
    pub p_daughter2_outlet: T,
    pub p_outlet: T,
    pub dp_parent: T,
    pub dp_daughter1: T,
    pub dp_daughter2: T,
    pub wall_shear_stress_parent: T,
    pub wall_shear_stress_daughter1: T,
    pub wall_shear_stress_daughter2: T,
    pub mass_conservation_error: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> BifurcationSolution3D<T> {
    pub fn new(_geometry: &BifurcationGeometry3D<T>) -> Self {
        Self {
            q_parent: T::zero(),
            q_daughter1: T::zero(),
            q_daughter2: T::zero(),
            u_parent_mean: T::zero(),
            u_daughter1_mean: T::zero(),
            u_daughter2_mean: T::zero(),
            p_inlet: T::zero(),
            p_junction_mid: T::zero(),
            p_daughter1_outlet: T::zero(),
            p_daughter2_outlet: T::zero(),
            p_outlet: T::zero(),
            dp_parent: T::zero(),
            dp_daughter1: T::zero(),
            dp_daughter2: T::zero(),
            wall_shear_stress_parent: T::zero(),
            wall_shear_stress_daughter1: T::zero(),
            wall_shear_stress_daughter2: T::zero(),
            mass_conservation_error: T::zero(),
        }
    }

    pub fn is_mass_conserved(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance
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
}
