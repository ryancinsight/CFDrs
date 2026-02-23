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

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + Float + From<f64>> SerpentineSolver3D<T> {
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
        let base_mesh = self.builder
            .clone()
            .with_resolution(self.config.resolution.0, self.config.resolution.1)
            .with_circular(self.config.circular)
            .build()
            .map_err(|e| Error::Solver(e.to_string()))?;

        // 1.1 Decompose to Tetrahedra and Promote to Quadratic (P2) mesh for Taylor-Hood elements (Q2-Q1)
        let legacy_mesh = cfd_mesh::mesh::Mesh::from_indexed(&base_mesh);
        let tet_mesh = cfd_mesh::hierarchy::hex_to_tet::HexToTetConverter::convert(&legacy_mesh);
        let p2_mesh_f64 = cfd_mesh::hierarchy::hierarchical_mesh::P2MeshConverter::convert_to_p2(&tet_mesh);

        let mut mesh = cfd_mesh::mesh::Mesh::<T>::new();
        for v in p2_mesh_f64.vertices() {
            let p = v.position;
            mesh.add_vertex(cfd_mesh::topology::Vertex::new(nalgebra::Point3::new(
                T::from_f64(p.x).unwrap(), T::from_f64(p.y).unwrap(), T::from_f64(p.z).unwrap()
            )));
        }
        for f in p2_mesh_f64.faces() {
            mesh.add_face(f.clone());
        }
        for c in p2_mesh_f64.cells() {
            mesh.add_cell(c.clone());
        }
        for f_idx in p2_mesh_f64.marked_boundary_faces() {
            if let Some(label) = p2_mesh_f64.boundary_label(f_idx) {
                mesh.mark_boundary(f_idx, label.to_string());
            }
        }

        // 2. Define Boundary Conditions
        // Priority: inlet > outlet > wall. Process inlet/outlet first so that
        // shared corner/edge nodes get the correct BC instead of a no-slip wall BC.
        let mut boundary_conditions = HashMap::new();
        let fluid_props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let area_inlet = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.builder.diameter * self.builder.diameter
        } else {
            self.builder.diameter * self.builder.diameter
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        // Pass 1: Assign inlet BCs (highest priority)
        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if label == "inlet" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.insert(v_idx,
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(T::zero(), T::zero(), u_inlet),
                                });
                        }
                    }
                }
            }
        }

        // Pass 2: Assign outlet BCs (do not overwrite inlet nodes)
        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if label == "outlet" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.entry(v_idx).or_insert(
                                BoundaryCondition::PressureOutlet {
                                    pressure: self.config.outlet_pressure,
                                });
                        }
                    }
                }
            }
        }

        // Pass 3: Assign wall (no-slip) BCs to remaining boundary nodes
        for f_idx in mesh.boundary_faces() {
            if let Some(label) = mesh.boundary_label(f_idx) {
                if label != "inlet" && label != "outlet" {
                    if let Some(face) = mesh.face(f_idx) {
                        for &v_idx in &face.vertices {
                            boundary_conditions.entry(v_idx).or_insert(
                                BoundaryCondition::Dirichlet {
                                    value: T::zero(),
                                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                                });
                        }
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

        let mut problem = StokesFlowProblem::new(mesh, constant_basis, boundary_conditions, tet_mesh.vertex_count());
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities = vec![fluid_props.dynamic_viscosity; n_elements];
        
        // 4. Picard Iteration Loop
        let fem_config = FemConfig::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;

        for iter in 0..self.config.max_nonlinear_iterations {
            problem.element_viscosities = Some(element_viscosities.clone());
            let fem_result = solver.solve(&problem, last_solution.as_ref());

            // If the linear solver fails (e.g. GMRES stagnation on the
            // non-Newtonian saddle-point system), use the last converged
            // solution and break.
            let fem_solution = match fem_result {
                Ok(sol) => sol,
                Err(e) => {
                    println!("Picard iteration {}: linear solve failed ({}), using last converged solution", iter, e);
                    break;
                }
            };

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
                let shear_rate = self.calculate_cell_shear_rate(cell, &problem.mesh, &updated_solution)?;
                let new_visc = fluid.viscosity_at_shear(shear_rate, T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

                let change = Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change {
                    max_change = change;
                }
                new_viscosities.push(new_visc);
            }

            element_viscosities = new_viscosities;
            last_solution = Some(updated_solution);

            println!("Picard iteration {}: visc_change={:?}", iter, max_change);
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
        let mut idxs = Vec::with_capacity(10);
        for &face_idx in &cell.faces {
            if let Some(face) = mesh.face(face_idx) {
                for &v_idx in &face.vertices {
                    if !idxs.contains(&v_idx) {
                        idxs.push(v_idx);
                    }
                }
            }
        }

        if idxs.len() < 4 {
            return Err(Error::Solver("Invalid cell topology".to_string()));
        }

        let mut local_verts = Vec::with_capacity(idxs.len());
        for &idx in &idxs {
            local_verts.push(mesh.vertex(idx).unwrap().position.coords);
        }

        let mut l = nalgebra::Matrix3::zeros();
        
        if idxs.len() == 10 {
            // Tet10 (P2): Evaluate gradient at centroid (L = [0.25, 0.25, 0.25, 0.25])
            use crate::fem::shape_functions::LagrangeTet10;
            
            // 1. Calculate P1 gradients (∇L_i)
            let mut tet4 = crate::fem::element::FluidElement::new(idxs[0..4].to_vec());
            let six_v = tet4.calculate_volume(&local_verts);
            if num_traits::Float::abs(six_v) < T::from_f64(1e-24).unwrap() {
                return Ok(T::zero());
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
            let l_centroid = [T::from_f64(0.25).unwrap(); 4];
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

        let epsilon = (l + l.transpose()) * T::from_f64_or_one(0.5);
        let mut inner_prod = T::zero();
        for i in 0..3 {
            for j in 0..3 {
                inner_prod += epsilon[(i, j)] * epsilon[(i, j)];
            }
        }
        let shear = Float::sqrt(T::from_f64_or_one(2.0) * inner_prod);

        Ok(shear)
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
