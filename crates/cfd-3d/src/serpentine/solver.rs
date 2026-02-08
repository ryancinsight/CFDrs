//! 3D FEM Navier-Stokes solver for Serpentine channels
//!
//! Solves the incompressible Navier-Stokes equations on 3D serpentine domains
//! using Finite Element Method with support for non-Newtonian blood rheology and
//! Dean flow analysis.

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::geometry::serpentine::SerpentineMeshBuilder;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D Serpentine solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineConfig3D<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
    for SerpentineConfig3D<T>
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
pub struct SerpentineSolver3D<T: RealField + Copy + Float> {
    builder: SerpentineMeshBuilder<T>,
    config: SerpentineConfig3D<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float> SerpentineSolver3D<T> {
    /// Create new solver from mesh builder and config
    pub fn new(builder: SerpentineMeshBuilder<T>, config: SerpentineConfig3D<T>) -> Self {
        Self { builder, config }
    }

    /// Solve Serpentine flow with given fluid
    pub fn solve<F: FluidTrait<T> + Clone>(
        &self,
        fluid: F,
    ) -> Result<SerpentineSolution3D<T>> {
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
            T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.builder.diameter * self.builder.diameter
        } else {
            self.builder.diameter * self.builder.diameter
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        boundary_conditions.entry(v_idx).or_insert_with(|| {
                            if label == "inlet" {
                                // For serpentine, the inlet might be slightly rotated.
                                // But StructuredGridBuilder used Z as axial, so it should be Z-aligned.
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
                                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero())]),
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
            let fem_solution = solver.solve(&problem).map_err(|e| Error::Solver(e.to_string()))?;
            
            let mut max_change = T::zero();
            let mut new_viscosities = Vec::with_capacity(n_elements);
            
            for (i, cell) in problem.mesh.cells().iter().enumerate() {
                let shear_rate = self.calculate_cell_shear_rate(cell, &problem.mesh, &fem_solution)?;
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
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;
        
        // 5. Extract Metrics
        let mut solution = SerpentineSolution3D::new();
        solution.u_inlet = u_inlet;
        solution.p_inlet = self.config.inlet_pressure;
        solution.p_outlet = self.config.outlet_pressure;
        solution.dp_total = solution.p_inlet - solution.p_outlet;
        
        // Calculate Dean Number: De = Re * sqrt(Dh / 2Rc)
        // For sine wave path x = A*sin(k*z), curvature kappa = |x''| / (1 + x'^2)^(3/2)
        // Max curvature at peaks: kappa_max = A*k^2. Radius Rc = 1/kappa_max = 1 / (A * (2pi/lambda)^2)
        let k = T::from_f64(2.0 * std::f64::consts::PI).unwrap() / self.builder.wavelength;
        let kappa_max = self.builder.amplitude * k * k;
        let rc = T::one() / Float::max(kappa_max, T::from_f64(1e-10).unwrap());
        
        let re = (fluid_props.density * u_inlet * self.builder.diameter) / fluid_props.dynamic_viscosity;
        solution.dean_number = re * Float::sqrt(self.builder.diameter / (T::from_f64(2.0).unwrap() * rc));

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

/// Complete solution to 3D Serpentine problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SerpentineSolution3D<T: RealField + Copy> {
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

impl<T: RealField + Copy> SerpentineSolution3D<T> {
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
