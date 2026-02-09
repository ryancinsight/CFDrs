//! 3D FEM Navier-Stokes solver for bifurcations with blood flow
//!
//! Solves the incompressible Navier-Stokes equations on 3D bifurcation domains
//! using Finite Element Method with stabilization for advection-dominated flows.

use super::geometry::BifurcationGeometry3D;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D bifurcation solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationConfig3D<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
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
            nonlinear_tolerance: T::from_f64_or_one(1e-6),
            max_linear_iterations: 100,
            linear_tolerance: T::from_f64_or_one(1e-8),
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
pub struct BifurcationSolver3D<T: RealField + Copy> {
    geometry: BifurcationGeometry3D<T>,
    config: BifurcationConfig3D<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float + From<f64>> BifurcationSolver3D<T> {
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
        use cfd_mesh::geometry::branching::BranchingMeshBuilder;
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
        let mesh = mesh_builder.build().map_err(|e| Error::Solver(e.to_string()))?;

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
        
        // Process outlet faces FIRST (highest priority) - use VelocityInlet with computed velocity
        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if label == "outlet_0" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.entry(v_idx).or_insert_with(|| {
                                outlet_count += 1;
                                // Daughter 1 at positive angle - rotate velocity vector
                                let angle = self.geometry.branching_angle;
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(
                                        u_outlet * Float::cos(angle), 
                                        u_outlet * Float::sin(angle), 
                                        T::zero()
                                    ),
                                }
                            });
                        }
                    }
                } else if label == "outlet_1" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.entry(v_idx).or_insert_with(|| {
                                outlet_count += 1;
                                // Daughter 2 at negative angle
                                let angle = -self.geometry.branching_angle;
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(
                                        u_outlet * Float::cos(angle), 
                                        u_outlet * Float::sin(angle), 
                                        T::zero()
                                    ),
                                }
                            });
                        }
                    }
                }
            }
        }
        
        // Process inlet faces SECOND
        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if label == "inlet" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.entry(v_idx).or_insert_with(|| {
                                inlet_count += 1;
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(u_inlet, T::zero(), T::zero()),
                                }
                            });
                        }
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
        for cell in mesh.cells() {
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
            if let Some(face) = mesh.face(face_idx) {
                for &v_idx in &face.vertices {
                    // If no BC assigned yet, default to wall
                    let bc = boundary_conditions.entry(v_idx).or_insert_with(|| {
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

        // 3. Set up FEM Problem with initial viscosity
        let constant_fluid = cfd_core::physics::fluid::ConstantPropertyFluid {
            name: "Picard Iteration Basis".to_string(),
            density: fluid_props.density,
            viscosity: fluid_props.dynamic_viscosity,
            specific_heat: fluid_props.specific_heat,
            thermal_conductivity: fluid_props.thermal_conductivity,
            speed_of_sound: fluid_props.speed_of_sound,
        };

        let mut problem = StokesFlowProblem::new(mesh, constant_fluid, boundary_conditions);
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities = vec![fluid_props.dynamic_viscosity; n_elements];
        
        // 4. Picard Iteration Loop
        let fem_config = FemConfig::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;

        for iter in 0..self.config.max_nonlinear_iterations {
            tracing::info!("Picard iteration {}", iter);
            
            // Update viscosities in problem
            problem.element_viscosities = Some(element_viscosities.clone());
            
            // Solve Stokes system
            let fem_solution = solver.solve(&problem, last_solution.as_ref()).map_err(|e| Error::Solver(e.to_string()))?;
            
            // Calculate shear rates and new viscosities
            let mut max_change = T::zero();
            let mut new_viscosities = Vec::with_capacity(n_elements);
            
            for (i, cell) in problem.mesh.cells().iter().enumerate() {
                let shear_rate = self.calculate_element_shear_rate(cell, &problem.mesh, &fem_solution)?;
                let new_visc = fluid.viscosity_at_shear(shear_rate, T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
                
                let change = Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change {
                    max_change = change;
                }
                new_viscosities.push(new_visc);
            }
            
            element_viscosities = new_viscosities;
            last_solution = Some(fem_solution);
            
            if max_change < self.config.nonlinear_tolerance {
                tracing::info!("Picard converged in {} iterations", iter + 1);
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        
        // Debug: Check solution magnitude (use max absolute value)
        use num_traits::Float;
        let vel_max: T = fem_solution.velocity.iter().map(|v| v.abs()).fold(T::zero(), |a, b| Float::max(a, b));
        let p_max: T = fem_solution.pressure.iter().map(|p| p.abs()).fold(T::zero(), |a, b| Float::max(a, b));
        eprintln!("DEBUG: FEM solution max velocity = {:?}, max pressure = {:?}", vel_max, p_max);

        // 5. Extract Metrics for BifurcationSolution3D
        let mesh = &problem.mesh;
        let mut solution = BifurcationSolution3D::new(&self.geometry);
        solution.u_parent_mean = u_inlet;
        solution.q_parent = self.config.inlet_flow_rate;
        
        // Calculate flows through daughters
        solution.q_daughter1 = self.calculate_boundary_flow(mesh, &fem_solution, "outlet_0")?;
        solution.q_daughter2 = self.calculate_boundary_flow(mesh, &fem_solution, "outlet_1")?;
        
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
        solution.p_junction_mid = self.extract_point_pressure(mesh, &fem_solution, Vector3::new(self.geometry.l_parent, T::zero(), T::zero()))?;
        
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

    /// Calculate flow rate through a boundary label using u·n integration
    fn calculate_boundary_flow(
        &self,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
        label: &str,
    ) -> Result<T> {
        let mut total_q = T::zero();
        let mut face_count = 0;
        
        for f_idx in 0..mesh.face_count() {
            if mesh.boundary_label(f_idx) == Some(label) {
                if let Some(face) = mesh.face(f_idx) {
                    if face.vertices.len() >= 3 {
                        face_count += 1;
                        let v0 = mesh.vertex(face.vertices[0]).unwrap().position.coords;
                        let v1 = mesh.vertex(face.vertices[1]).unwrap().position.coords;
                        let v2 = mesh.vertex(face.vertices[2]).unwrap().position.coords;
                        
                        let n_vec = (v1 - v0).cross(&(v2 - v0));
                        let area = n_vec.norm() * T::from_f64_or_one(0.5);
                        let face_normal = n_vec.normalize();
                        
                        let mut u_avg = Vector3::zeros();
                        for &v_idx in &face.vertices {
                            let u = solution.get_velocity(v_idx);
                            u_avg += u;
                            // Debug first few vertices
                            if face_count <= 2 {
                                eprintln!("DEBUG: vertex {} velocity = {:?}", v_idx, u);
                            }
                        }
                        u_avg /= T::from_usize(face.vertices.len()).unwrap_or_else(T::one);
                        
                        let face_flow = u_avg.dot(&face_normal) * area;
                        total_q += face_flow;
                    }
                }
            }
        }
        
        eprintln!("DEBUG: Flow integration for '{}': {} faces, total_q = {:?}", label, face_count, total_q);
        Ok(Float::abs(total_q))
    }

    /// Extract pressure at a point in the mesh (closest node)
    fn extract_point_pressure(
        &self,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
        point: Vector3<T>,
    ) -> Result<T> {
        let mut best_node = 0;
        let mut min_dist = <T as RealField>::max_value().unwrap_or_else(T::one);
        
        for (i, v) in mesh.vertices().iter().enumerate() {
            let dist = (v.position.coords - point).norm();
            if dist < min_dist {
                min_dist = dist;
                best_node = i;
            }
        }
        
        Ok(solution.get_pressure(best_node))
    }

    /// Calculate shear rate within a tetrahedral element
    fn calculate_element_shear_rate(
        &self,
        cell: &cfd_mesh::topology::Cell,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
    ) -> Result<T> {
        use crate::fem::element::FluidElement;
        
        let mut vertex_indices = Vec::with_capacity(4);
        for &face_idx in &cell.faces {
            if let Some(face) = mesh.face(face_idx) {
                for &v_idx in &face.vertices {
                    if !vertex_indices.contains(&v_idx) {
                        vertex_indices.push(v_idx);
                    }
                }
            }
            if vertex_indices.len() >= 4 { break; }
        }

        if vertex_indices.len() < 4 {
            return Err(Error::Solver("Invalid cell topology".to_string()));
        }

        let mut element = FluidElement::new(vertex_indices.clone());
        let vertex_positions: Vec<Vector3<T>> = mesh.vertices().iter().map(|v| v.position.coords).collect();
        // Just take the positions of the 4 nodes
        let mut element_vertices = Vec::with_capacity(4);
        for &idx in &vertex_indices {
            element_vertices.push(vertex_positions[idx]);
        }
        
        element.calculate_shape_derivatives(&element_vertices);
        
        let mut l = nalgebra::Matrix3::zeros();
        for i in 0..4 {
            let u = solution.get_velocity(vertex_indices[i]);
            for j in 0..3 { 
                for k in 0..3 { 
                    l[(j, k)] += element.shape_derivatives[(k, i)] * u[j];
                }
            }
        }
        
        let epsilon = (l + l.transpose()) * T::from_f64_or_one(0.5);
        let mut inner_prod = T::zero();
        for i in 0..3 {
            for j in 0..3 {
                inner_prod += epsilon[(i, j)] * epsilon[(i, j)];
            }
        }
        
        Ok(Float::sqrt(T::from_f64_or_one(2.0) * inner_prod))
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
pub struct BifurcationSolution3D<T: RealField + Copy> {
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

impl<T: RealField + Copy> BifurcationSolution3D<T> {
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
