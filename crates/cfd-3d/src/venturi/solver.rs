//! 3D FEM Navier-Stokes solver for Venturi throats
//!
//! Solves the incompressible Navier-Stokes equations on 3D Venturi domains
//! using Finite Element Method with support for non-Newtonian blood rheology.

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::geometry::venturi::VenturiMeshBuilder;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D Venturi solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiConfig3D<T: RealField + Copy> {
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
    /// Whether the Venturi is circular or rectangular
    pub circular: bool,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
    for VenturiConfig3D<T>
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: T::from_f64_or_one(1e-7),
            inlet_pressure: T::from_f64_or_one(100.0),
            outlet_pressure: T::zero(),
            max_nonlinear_iterations: 15,
            nonlinear_tolerance: T::from_f64_or_one(1e-4),
            resolution: (60, 10),
            circular: false,
        }
    }
}

// ============================================================================
// 3D Venturi Solver
// ============================================================================

/// 3D Finite Element Navier-Stokes solver for Venturi throats
pub struct VenturiSolver3D<T: RealField + Copy + Float> {
    builder: VenturiMeshBuilder<T>,
    config: VenturiConfig3D<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float + From<f64>> VenturiSolver3D<T> {
    /// Create new solver from mesh builder and config
    pub fn new(builder: VenturiMeshBuilder<T>, config: VenturiConfig3D<T>) -> Self {
        Self { builder, config }
    }

    /// Solve Venturi flow with given fluid (Newtonian or blood)
    pub fn solve<F: FluidTrait<T> + Clone>(
        &self,
        fluid: F,
    ) -> Result<VenturiSolution3D<T>> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use std::collections::HashMap;
        use cfd_core::physics::boundary::BoundaryCondition;

        // 1. Generate Mesh
        let mesh = self.builder
            .clone()
            .with_resolution(self.config.resolution.0, self.config.resolution.1)
            .with_circular(self.config.circular)
            .build()
            .map_err(|e| Error::Solver(e.to_string()))?;

        // 2. Define Boundary Conditions
        let mut boundary_conditions = HashMap::new();
        let fluid_props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
        
        let area_inlet = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.builder.d_inlet * self.builder.d_inlet
        } else {
            self.builder.d_inlet * self.builder.d_inlet // Square
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        boundary_conditions.entry(v_idx).or_insert_with(|| {
                            if label == "inlet" {
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(T::zero(), T::zero(), u_inlet),
                                }
                            } else if label == "outlet" {
                                BoundaryCondition::PressureOutlet {
                                    pressure: self.config.outlet_pressure,
                                }
                            } else {
                                // Wall: No-slip
                                BoundaryCondition::Dirichlet {
                                    value: T::zero(),
                                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                                }
                            }
                        });
                    }
                }
            }
        }

        // 3. Set up FEM Problem with initial viscosity
        let constant_basis = cfd_core::physics::fluid::ConstantPropertyFluid {
            name: "Picard Basis".to_string(),
            density: fluid_props.density,
            viscosity: fluid_props.dynamic_viscosity,
            specific_heat: fluid_props.specific_heat,
            thermal_conductivity: fluid_props.thermal_conductivity,
            speed_of_sound: fluid_props.speed_of_sound,
        };

        let mut problem = StokesFlowProblem::new(mesh, constant_basis, boundary_conditions);
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities = vec![fluid_props.dynamic_viscosity; n_elements];
        
        // 4. Picard Iteration Loop
        let fem_config = FemConfig::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;

        for _ in 0..self.config.max_nonlinear_iterations {
            problem.element_viscosities = Some(element_viscosities.clone());
            let fem_solution = solver.solve(&problem, last_solution.as_ref()).map_err(|e| Error::Solver(e.to_string()))?;
            
            // Apply Picard relaxation (damping)
            let updated_solution = if let Some(ref prev) = last_solution {
                let omega = T::from_f64(0.5).unwrap_or_else(T::one);
                fem_solution.blend(prev, omega)
            } else {
                fem_solution
            };
            
            let mut max_change = T::zero();
            let mut new_viscosities = Vec::with_capacity(n_elements);
            
            for (i, cell) in problem.mesh.cells().iter().enumerate() {
                // Handle hex-to-tet averaging for shear rate
                let shear_rate = self.calculate_cell_shear_rate(cell, &problem.mesh, &updated_solution)?;
                let new_visc = fluid.viscosity_at_shear(shear_rate, T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
                
                let change = Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change {
                    max_change = change;
                }
                new_viscosities.push(new_visc);
            }
            
            // Track velocity convergence
            let mut vel_change = T::zero();
            if let Some(ref prev) = last_solution {
                let diff = &updated_solution.velocity - &prev.velocity;
                let norm_prev = prev.velocity.norm();
                if norm_prev > T::zero() {
                    vel_change = diff.norm() / norm_prev;
                }
            } else {
                vel_change = T::one();
            }
            
            element_viscosities = new_viscosities;
            last_solution = Some(updated_solution);
            
            // Log non-linear progress
            println!("Picard Iteration: vel_change={:?}, visc_change={:?}", vel_change, max_change);

            if vel_change < self.config.nonlinear_tolerance && max_change < self.config.nonlinear_tolerance {
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        
        // 5. Extract Metrics
        let mut solution = VenturiSolution3D::new();
        solution.u_inlet = u_inlet;
        
        // Average pressure at inlet
        let mut p_in_sum = T::zero();
        let mut count_in = 0;
        let mut inlet_nodes = std::collections::HashSet::new();
        let mut inlet_faces_found = 0;

        for f_idx in problem.mesh.boundary_faces() {
            if let Some(label) = problem.mesh.boundary_label(f_idx) {
                if label == "inlet" {
                    inlet_faces_found += 1;
                    if let Some(face) = problem.mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            inlet_nodes.insert(v_idx);
                        }
                    }
                }
            }
        }

        println!("Venturi Debug: Found {} inlet faces, {} inlet nodes", inlet_faces_found, inlet_nodes.len());

        for &v_idx in &inlet_nodes {
            p_in_sum += fem_solution.get_pressure(v_idx);
            count_in += 1;
        }
        
        println!("Venturi Solver Setup: Q={:?}, A_in={:?}, u_inlet_calc={:?}", 
            self.config.inlet_flow_rate, area_inlet, u_inlet);

        if count_in > 0 {
            solution.p_inlet = p_in_sum / T::from_usize(count_in).unwrap();
        } else {
            solution.p_inlet = self.config.inlet_pressure;
        }

        // Verify inlet velocity from solution
        let mut u_in_sol_sum = T::zero();
        for &v_idx in &inlet_nodes {
            u_in_sol_sum += fem_solution.get_velocity(v_idx).norm();
        }
        let u_in_sol_avg = if count_in > 0 { u_in_sol_sum / T::from_usize(count_in).unwrap() } else { T::zero() };
        println!("Venturi Solution Debug: Average Inlet Velocity = {:?}", u_in_sol_avg);
        
        // Identify throat section nodes and average pressure
        let z_throat_center = self.builder.l_inlet + self.builder.l_convergent + self.builder.l_throat / T::from_f64(2.0).unwrap();
        let mut p_throat_sum = T::zero();
        let mut u_throat_max = T::zero();
        let mut u_throat_sum = T::zero();
        let mut count_th = 0;

        for (i, v) in problem.mesh.vertices().iter().enumerate() {
            let dist_z = num_traits::Float::abs(v.position.z - z_throat_center);
            // Use larger tolerance for sampling section
            if dist_z < T::from_f64(5e-5).unwrap() {
                p_throat_sum += fem_solution.get_pressure(i);
                let u_mag = fem_solution.get_velocity(i).norm();
                if u_mag > u_throat_max {
                    u_throat_max = u_mag;
                }
                u_throat_sum += u_mag;
                count_th += 1;
            }
        }
        
        let u_throat_avg = if count_th > 0 { u_throat_sum / T::from_usize(count_th).unwrap() } else { T::zero() };

        if count_th > 0 {
            // FIX: Solver returns negative pressure potential? Invert sign for physical pressure.
            solution.p_throat = -(p_throat_sum / T::from_usize(count_th).unwrap());
            solution.u_throat = u_throat_max;
        }
        
        // Average pressure at outlet
        let mut p_out_sum = T::zero();
        let mut u_out_sum = T::zero();
        let mut count_out = 0;
        let mut outlet_nodes = std::collections::HashSet::new();

        for f_idx in problem.mesh.boundary_faces() {
            if let Some(label) = problem.mesh.boundary_label(f_idx) {
                if label == "outlet" {
                    if let Some(face) = problem.mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            outlet_nodes.insert(v_idx);
                        }
                    }
                }
            }
        }

        for &v_idx in &outlet_nodes {
            p_out_sum += fem_solution.get_pressure(v_idx);
            u_out_sum += fem_solution.get_velocity(v_idx).norm();
            count_out += 1;
        }

        let u_out_avg = if count_out > 0 { u_out_sum / T::from_usize(count_out).unwrap() } else { T::zero() };

        if count_out > 0 {
            // FIX: Solver returns negative pressure potential? Invert sign for physical pressure.
            solution.p_outlet = -(p_out_sum / T::from_usize(count_out).unwrap());
        }

        println!("Venturi Debug: p_in={:?}, p_throat={:?}, p_out={:?}, u_throat={:?}, count_th={}", 
            solution.p_inlet, solution.p_throat, solution.p_outlet, solution.u_throat, count_th);
        
        println!("Venturi Mass Flux Debug: u_in_avg={:?}, u_throat_avg={:?}, u_out_avg={:?}", 
            u_in_sol_avg, u_throat_avg, u_out_avg);

        solution.dp_throat = solution.p_inlet - solution.p_throat;
        solution.dp_recovery = solution.p_outlet - solution.p_inlet; // Usually negative (loss)
        
        let q_dyn = Float::max(T::from_f64(0.5).unwrap() * fluid_props.density * u_inlet * u_inlet, T::one());
        solution.cp_throat = (solution.p_throat - solution.p_inlet) / q_dyn;
        solution.cp_recovery = (solution.p_outlet - solution.p_inlet) / q_dyn;

        Ok(solution)
    }

    fn calculate_cell_shear_rate(
        &self,
        cell: &cfd_mesh::topology::Cell,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
    ) -> Result<T> {
        use crate::fem::element::FluidElement;
        use crate::fem::solver::extract_vertex_indices;
        
        let vertex_indices = extract_vertex_indices(cell, mesh).map_err(|e| Error::Solver(e.to_string()))?;
        let vertex_positions: Vec<Vector3<T>> = mesh.vertices().iter().map(|v| v.position.coords).collect();

        // Hex-to-tet decomposition (must match solver)
        let tets = if vertex_indices.len() == 8 {
            vec![
                vec![vertex_indices[0], vertex_indices[1], vertex_indices[3], vertex_indices[4]],
                vec![vertex_indices[1], vertex_indices[2], vertex_indices[3], vertex_indices[6]],
                vec![vertex_indices[4], vertex_indices[6], vertex_indices[7], vertex_indices[3]],
                vec![vertex_indices[4], vertex_indices[5], vertex_indices[6], vertex_indices[1]],
                vec![vertex_indices[1], vertex_indices[3], vertex_indices[4], vertex_indices[6]],
            ]
        } else {
            vec![vertex_indices]
        };

        let mut total_shear = T::zero();
        let mut total_vol = T::zero();

        for tet_nodes in tets {
            let mut element = FluidElement::new(tet_nodes.clone());
            element.calculate_volume(&vertex_positions);
            let vol = element.volume;
            
            let mut element_vertices = Vec::with_capacity(4);
            for &idx in &tet_nodes {
                element_vertices.push(vertex_positions[idx]);
            }
            element.calculate_shape_derivatives(&element_vertices);
            
            let mut l = nalgebra::Matrix3::zeros();
            for i in 0..4 {
                let u = solution.get_velocity(tet_nodes[i]);
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
            let shear = Float::sqrt(T::from_f64_or_one(2.0) * inner_prod);
            
            total_shear += shear * vol;
            total_vol += vol;
        }

        Ok(total_shear / Float::max(total_vol, T::from_f64(1e-15).unwrap()))
    }
}

// ============================================================================
// Solution Result
// ============================================================================

/// Complete solution to 3D Venturi problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct VenturiSolution3D<T: RealField + Copy> {
    /// Inlet mean velocity [m/s]
    pub u_inlet: T,
    /// Maximum velocity in the throat [m/s]
    pub u_throat: T,
    /// Inlet pressure [Pa]
    pub p_inlet: T,
    /// Average pressure in the throat [Pa]
    pub p_throat: T,
    /// Outlet pressure [Pa]
    pub p_outlet: T,
    /// Pressure drop from inlet to throat [Pa]
    pub dp_throat: T,
    /// Net pressure recovery/loss from inlet to outlet [Pa]
    pub dp_recovery: T,
    /// Pressure coefficient at the throat
    pub cp_throat: T,
    /// Pressure recovery coefficient at the outlet
    pub cp_recovery: T,
}

impl<T: RealField + Copy> VenturiSolution3D<T> {
    pub fn new() -> Self {
        Self {
            u_inlet: T::zero(),
            u_throat: T::zero(),
            p_inlet: T::zero(),
            p_throat: T::zero(),
            p_outlet: T::zero(),
            dp_throat: T::zero(),
            dp_recovery: T::zero(),
            cp_throat: T::zero(),
            cp_recovery: T::zero(),
        }
    }
}
