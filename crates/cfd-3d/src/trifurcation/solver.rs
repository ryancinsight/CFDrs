//! 3D solver for flow in trifurcations
//!
//! Provides simulation for three-way branching vessels with non-Newtonian blood.

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::{Fluid as FluidTrait, NonNewtonianFluid};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};
use crate::trifurcation::geometry::TrifurcationGeometry3D;

/// Configuration for 3D trifurcation solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationConfig3D<T: RealField + Copy> {
    pub inlet_flow_rate: T,
    pub inlet_pressure: T,
    pub outlet_pressures: [T; 3],
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default for TrifurcationConfig3D<T> {
    fn default() -> Self {
        Self {
            inlet_flow_rate: T::from_f64_or_one(1e-8),
            inlet_pressure: T::from_f64_or_one(100.0),
            outlet_pressures: [T::zero(), T::zero(), T::zero()],
        }
    }
}

/// Solution result for 3D trifurcation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationSolution3D<T: RealField + Copy> {
    pub flow_rates: [T; 4], // [parent, d1, d2, d3]
    pub mean_velocities: [T; 4],
    pub wall_shear_stresses: [T; 4],
    pub pressure_drops: [T; 4],
    pub mass_conservation_error: T,
}

/// 3D Trifurcation flow solver
pub struct TrifurcationSolver3D<T: RealField + Copy> {
    pub geometry: TrifurcationGeometry3D<T>,
    pub config: TrifurcationConfig3D<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float> TrifurcationSolver3D<T> {
    pub fn new(geometry: TrifurcationGeometry3D<T>, config: TrifurcationConfig3D<T>) -> Self {
        Self { geometry, config }
    }

    /// Solve trifurcation flow
    pub fn solve<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        fluid: F,
    ) -> Result<TrifurcationSolution3D<T>> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use cfd_mesh::geometry::branching::BranchingMeshBuilder;
        use std::collections::HashMap;
        use cfd_core::physics::boundary::BoundaryCondition;

        // 1. Generate Mesh
        let mesh_builder = BranchingMeshBuilder::trifurcation(
            self.geometry.d_parent,
            self.geometry.l_parent,
            self.geometry.d_daughters[0],
            self.geometry.l_daughters[0],
            self.geometry.branching_angles[0],
            8, // resolution factor
        );
        let mesh = mesh_builder.build().map_err(|e| Error::Solver(e.to_string()))?;

        // 2. Define Boundary Conditions
        let mut boundary_conditions = HashMap::new();
        let fluid_props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
        
        let inlet_area = T::from_f64_or_one(std::f64::consts::PI / 4.0) * num_traits::Float::powf(self.geometry.d_parent, T::from_f64_or_one(2.0));
        let u_inlet = self.config.inlet_flow_rate / inlet_area;

        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        let bc = match label {
                            "inlet" => BoundaryCondition::Dirichlet {
                                value: u_inlet,
                                component_values: Some(vec![Some(u_inlet), Some(T::zero()), Some(T::zero())]),
                            },
                            "outlet_0" | "outlet_1" | "outlet_2" => BoundaryCondition::Neumann {
                                gradient: self.config.outlet_pressures[0],
                            },
                            "wall" => BoundaryCondition::Dirichlet {
                                value: T::zero(),
                                component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero())]),
                            },
                            _ => continue,
                        };
                        boundary_conditions.insert(v_idx, bc);
                    }
                }
            }
        }

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

        for iter in 0..20 { 
            problem.element_viscosities = Some(element_viscosities.clone());
            let fem_solution = solver.solve(&problem).map_err(|e| Error::Solver(e.to_string()))?;
            
            let mut max_change = T::zero();
            let mut new_viscosities = Vec::with_capacity(n_elements);
            
            for (i, cell) in problem.mesh.cells().iter().enumerate() {
                let shear_rate = self.calculate_element_shear_rate(cell, &problem.mesh, &fem_solution)?;
                let new_visc = fluid.apparent_viscosity(shear_rate);
                
                let change = Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change {
                    max_change = change;
                }
                new_viscosities.push(new_visc);
            }
            
            element_viscosities = new_viscosities;
            last_solution = Some(fem_solution);
            
            if max_change < T::from_f64_or_one(1e-6) {
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        let mesh = &problem.mesh;

        // 5. Build Result
        let q_parent = self.config.inlet_flow_rate;
        let q_d1 = self.calculate_boundary_flow(mesh, &fem_solution, "outlet_0")?;
        let q_d2 = self.calculate_boundary_flow(mesh, &fem_solution, "outlet_1")?;
        let q_d3 = self.calculate_boundary_flow(mesh, &fem_solution, "outlet_2")?;
        
        let mut solution = TrifurcationSolution3D {
            flow_rates: [q_parent, q_d1, q_d2, q_d3],
            mean_velocities: [u_inlet, T::zero(), T::zero(), T::zero()], 
            wall_shear_stresses: [T::zero(), T::zero(), T::zero(), T::zero()],
            pressure_drops: [T::zero(), T::zero(), T::zero(), T::zero()],
            mass_conservation_error: Float::abs(q_parent - (q_d1 + q_d2 + q_d3)),
        };

        // Fill in mean velocities for daughters
        let a_d = T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.geometry.d_daughters[0] * self.geometry.d_daughters[0];
        solution.mean_velocities[1] = q_d1 / a_d;
        solution.mean_velocities[2] = q_d2 / a_d;
        solution.mean_velocities[3] = q_d3 / a_d;

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
        
        for f_idx in 0..mesh.face_count() {
            if mesh.boundary_label(f_idx) == Some(label) {
                if let Some(face) = mesh.face(f_idx) {
                    if face.vertices.len() >= 3 {
                        let v0 = mesh.vertex(face.vertices[0]).unwrap().position.coords;
                        let v1 = mesh.vertex(face.vertices[1]).unwrap().position.coords;
                        let v2 = mesh.vertex(face.vertices[2]).unwrap().position.coords;
                        
                        let n_vec = (v1 - v0).cross(&(v2 - v0));
                        let area = n_vec.norm() * T::from_f64_or_one(0.5);
                        let face_normal = n_vec.normalize();
                        
                        let mut u_avg = Vector3::zeros();
                        for &v_idx in &face.vertices {
                            u_avg += solution.get_velocity(v_idx);
                        }
                        u_avg /= T::from_usize(face.vertices.len()).unwrap_or_else(T::one);
                        
                        total_q += u_avg.dot(&face_normal) * area;
                    }
                }
            }
        }
        
        Ok(Float::abs(total_q))
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
}
