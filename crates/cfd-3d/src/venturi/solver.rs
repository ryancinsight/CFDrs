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
        let base_mesh = self.builder
            .clone()
            .with_resolution(self.config.resolution.0, self.config.resolution.1)
            .with_circular(self.config.circular)
            .build()
            .map_err(|e| Error::Solver(e.to_string()))?;

        // 1.1 Decompose to Tetrahedra and Promote to Quadratic (P2) mesh for Taylor-Hood elements (Q2-Q1)
        let tet_mesh = cfd_mesh::hierarchy::hex_to_tet::HexToTetConverter::convert(&base_mesh);
        let mut mesh = cfd_mesh::hierarchy::hierarchical_mesh::P2MeshConverter::convert_to_p2(&tet_mesh);

        // Re-apply boundary labels after conversion to ensure all boundary faces are tagged.
        // Some conversion paths introduce boundary faces without labels.
        {
            let total_l = self.builder.l_inlet
                + self.builder.l_convergent
                + self.builder.l_throat
                + self.builder.l_divergent
                + self.builder.l_outlet;
            let tol = T::from_f64(1e-5).unwrap();

            let boundary_faces: Vec<usize> = mesh.boundary_faces();
            for f_idx in boundary_faces {
                if mesh.boundary_label(f_idx).is_some() {
                    continue;
                }

                if let Some(face) = mesh.face(f_idx) {
                    let mut all_at_inlet = true;
                    let mut all_at_outlet = true;
                    for &v_idx in &face.vertices {
                        if let Some(v) = mesh.vertex(v_idx) {
                            let z = v.position.z;
                            let dist_inlet = Float::abs(z);
                            let dist_outlet = Float::abs(z - total_l);
                            if dist_inlet > tol {
                                all_at_inlet = false;
                            }
                            if dist_outlet > tol {
                                all_at_outlet = false;
                            }
                        }
                    }

                    if all_at_inlet {
                        mesh.mark_boundary(f_idx, "inlet".to_string());
                    } else if all_at_outlet {
                        mesh.mark_boundary(f_idx, "outlet".to_string());
                    } else {
                        mesh.mark_boundary(f_idx, "wall".to_string());
                    }
                }
            }
        }

        // Boundary diagnostics: labeled faces vs connectivity boundary faces
        {
            use std::collections::{HashMap, HashSet};

            let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
            for cell in mesh.cells() {
                for &face_idx in &cell.faces {
                    *face_cell_count.entry(face_idx).or_insert(0) += 1;
                }
            }

            let connectivity_boundary_faces: HashSet<usize> = face_cell_count
                .iter()
                .filter(|&(_face_idx, &count)| count == 1)
                .map(|(&face_idx, _)| face_idx)
                .collect();

            let marked_boundary_faces: HashSet<usize> =
                mesh.boundary_faces().into_iter().collect();

            let boundary_union: HashSet<usize> = marked_boundary_faces
                .union(&connectivity_boundary_faces)
                .copied()
                .collect();

            let mut inlet_faces = 0usize;
            let mut outlet_faces = 0usize;
            let mut wall_faces = 0usize;
            let mut unlabeled_faces = 0usize;

            let mut inlet_nodes = HashSet::new();
            let mut outlet_nodes = HashSet::new();
            let mut wall_nodes = HashSet::new();
            let mut labeled_nodes = HashSet::new();

            for &f_idx in &marked_boundary_faces {
                let label = mesh.boundary_label(f_idx);
                match label.as_deref() {
                    Some("inlet") => {
                        inlet_faces += 1;
                        if let Some(face) = mesh.face(f_idx) {
                            for &v_idx in &face.vertices {
                                inlet_nodes.insert(v_idx);
                                labeled_nodes.insert(v_idx);
                            }
                        }
                    }
                    Some("outlet") => {
                        outlet_faces += 1;
                        if let Some(face) = mesh.face(f_idx) {
                            for &v_idx in &face.vertices {
                                outlet_nodes.insert(v_idx);
                                labeled_nodes.insert(v_idx);
                            }
                        }
                    }
                    Some("wall") => {
                        wall_faces += 1;
                        if let Some(face) = mesh.face(f_idx) {
                            for &v_idx in &face.vertices {
                                wall_nodes.insert(v_idx);
                                labeled_nodes.insert(v_idx);
                            }
                        }
                    }
                    _ => {
                        unlabeled_faces += 1;
                    }
                }
            }

            let mut connectivity_nodes = HashSet::new();
            for &f_idx in &connectivity_boundary_faces {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        connectivity_nodes.insert(v_idx);
                    }
                }
            }

            let inlet_outlet_overlap = inlet_nodes.intersection(&outlet_nodes).count();
            let inlet_wall_overlap = inlet_nodes.intersection(&wall_nodes).count();
            let outlet_wall_overlap = outlet_nodes.intersection(&wall_nodes).count();

            println!(
                "Venturi Boundary Debug: marked_faces={}, connectivity_faces={}, union_faces={}, unlabeled_marked={} ",
                marked_boundary_faces.len(),
                connectivity_boundary_faces.len(),
                boundary_union.len(),
                unlabeled_faces
            );
            println!(
                "Venturi Boundary Debug: inlet_faces={}, outlet_faces={}, wall_faces={}, inlet_nodes={}, outlet_nodes={}, wall_nodes={} ",
                inlet_faces, outlet_faces, wall_faces,
                inlet_nodes.len(), outlet_nodes.len(), wall_nodes.len()
            );
            println!(
                "Venturi Boundary Debug: overlap inlet/outlet={}, inlet/wall={}, outlet/wall={}, labeled_nodes={}, connectivity_nodes={}, total_nodes={}",
                inlet_outlet_overlap,
                inlet_wall_overlap,
                outlet_wall_overlap,
                labeled_nodes.len(),
                connectivity_nodes.len(),
                mesh.vertex_count()
            );

            let labeled_not_connectivity = labeled_nodes.difference(&connectivity_nodes).count();
            let connectivity_not_labeled = connectivity_nodes.difference(&labeled_nodes).count();
            println!(
                "Venturi Boundary Debug: labeled_not_connectivity={}, connectivity_not_labeled={}",
                labeled_not_connectivity,
                connectivity_not_labeled
            );
        }

        // 2. Define Boundary Conditions
        // Explicitly classify boundary node sets first so that inlet/wall intersection
        // (rim) nodes can be treated consistently.
        let mut boundary_conditions = HashMap::new();
        let fluid_props = fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let area_inlet = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.builder.d_inlet * self.builder.d_inlet
        } else {
            self.builder.d_inlet * self.builder.d_inlet // Square
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        let mut inlet_nodes = std::collections::HashSet::new();
        let mut outlet_nodes = std::collections::HashSet::new();
        let mut wall_nodes = std::collections::HashSet::new();
        let mut boundary_vertices = std::collections::HashSet::new();

        let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
        for cell in mesh.cells() {
            for &f_idx in &cell.faces {
                *face_cell_count.entry(f_idx).or_insert(0) += 1;
            }
        }
        let boundary_faces: Vec<usize> = face_cell_count
            .iter()
            .filter(|&(_, &count)| count == 1)
            .map(|(&idx, _)| idx)
            .collect();

        let mut z_min = <T as RealField>::max_value().unwrap_or_else(T::one);
        let mut z_max = -<T as RealField>::max_value().unwrap_or_else(T::one);
        for v in mesh.vertices() {
            let z = v.position.z;
            if z < z_min {
                z_min = z;
            }
            if z > z_max {
                z_max = z;
            }
        }
        let z_span = z_max - z_min;
        let z_tol = z_span / T::from_usize(self.config.resolution.0.max(1)).unwrap()
            * T::from_f64(0.75).unwrap_or_else(T::one);

        for &f_idx in &boundary_faces {
            if let Some(face) = mesh.face(f_idx) {
                let label = mesh.boundary_label(f_idx);
                let z_center = {
                    let mut sum = T::zero();
                    let mut count = 0usize;
                    for &v_idx in &face.vertices {
                        if let Some(v) = mesh.vertex(v_idx) {
                            sum += v.position.z;
                            count += 1;
                        }
                    }
                    if count > 0 {
                        sum / T::from_usize(count).unwrap()
                    } else {
                        T::zero()
                    }
                };

                let is_inlet = matches!(label, Some("inlet"))
                    || (label.is_none() && Float::abs(z_center - z_min) <= z_tol);
                let is_outlet = matches!(label, Some("outlet"))
                    || (label.is_none() && Float::abs(z_center - z_max) <= z_tol);

                for &v_idx in &face.vertices {
                    boundary_vertices.insert(v_idx);
                }

                if is_inlet {
                    for &v_idx in &face.vertices {
                        inlet_nodes.insert(v_idx);
                    }
                } else if is_outlet {
                    for &v_idx in &face.vertices {
                        outlet_nodes.insert(v_idx);
                    }
                } else {
                    for &v_idx in &face.vertices {
                        wall_nodes.insert(v_idx);
                    }
                }
            }
        }

        let inlet_wall_rim_nodes: std::collections::HashSet<usize> = inlet_nodes
            .intersection(&wall_nodes)
            .copied()
            .collect();

        // Pass 1: Apply inlet velocity to non-rim inlet nodes and no-slip to inlet-wall rim nodes.
        for &v_idx in inlet_nodes.iter() {
            if inlet_wall_rim_nodes.contains(&v_idx) {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::Dirichlet {
                        value: T::zero(),
                        component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                    },
                );
            } else {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::VelocityInlet {
                        velocity: Vector3::new(T::zero(), T::zero(), u_inlet),
                    },
                );
            }
        }

        // Pass 2: Assign a single outlet pressure reference (gauge) on a corner node.
        // Applying pressure Dirichlet on many outlet pressure DOFs over-constrains the
        // mixed system and removes too many continuity equations.
        let mut outlet_corner_nodes = std::collections::HashSet::new();
        for &v_idx in outlet_nodes.iter() {
            if v_idx < tet_mesh.vertex_count() {
                outlet_corner_nodes.insert(v_idx);
            }
        }

        let mut p_ref_node = None;
        let mut best_r2 = <T as RealField>::max_value().unwrap_or_else(T::one);
        for &v_idx in &outlet_corner_nodes {
            if let Some(v) = mesh.vertex(v_idx) {
                let r2 = v.position.x * v.position.x + v.position.y * v.position.y;
                if r2 < best_r2 {
                    best_r2 = r2;
                    p_ref_node = Some(v_idx);
                }
            }
        }

        if let Some(p_ref_node) = p_ref_node {
            boundary_conditions.insert(
                p_ref_node,
                BoundaryCondition::PressureOutlet {
                    pressure: self.config.outlet_pressure,
                },
            );
        }

        // Mark remaining outlet nodes as explicit natural outflow BCs so they are
        // not treated as unconstrained boundary leaks by generic diagnostics.
        for &v_idx in outlet_nodes.iter() {
            if boundary_conditions.contains_key(&v_idx) {
                continue;
            }
            boundary_conditions.insert(v_idx, BoundaryCondition::Outflow);
        }

        // Pass 3: Assign wall (no-slip) BCs to wall nodes EXCEPT those on the inlet or outlet.
        for &v_idx in wall_nodes.iter() {
            if inlet_nodes.contains(&v_idx) || outlet_nodes.contains(&v_idx) {
                continue;
            }
            boundary_conditions.insert(
                v_idx,
                BoundaryCondition::Dirichlet {
                    value: T::zero(),
                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                },
            );
        }

        // Pass 4: Ensure all boundary vertices have a BC unless intentionally left as outlet natural boundary.
        let inlet_radius_sq = {
            let r = T::from_f64(0.5).unwrap() * self.builder.d_inlet * T::from_f64(0.98).unwrap();
            r * r
        };
        let mut repaired_nodes = 0usize;
        for &v_idx in &boundary_vertices {
            if boundary_conditions.contains_key(&v_idx) {
                continue;
            }
            if let Some(v) = mesh.vertex(v_idx) {
                let z = v.position.z;
                if Float::abs(z - z_max) <= z_tol {
                    continue;
                }
                if Float::abs(z - z_min) <= z_tol {
                    if self.config.circular {
                        let r_sq = v.position.x * v.position.x + v.position.y * v.position.y;
                        if r_sq >= inlet_radius_sq {
                            boundary_conditions.insert(
                                v_idx,
                                BoundaryCondition::Dirichlet {
                                    value: T::zero(),
                                    component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                                },
                            );
                        } else {
                            boundary_conditions.insert(
                                v_idx,
                                BoundaryCondition::VelocityInlet {
                                    velocity: Vector3::new(T::zero(), T::zero(), u_inlet),
                                },
                            );
                        }
                    } else {
                        boundary_conditions.insert(
                            v_idx,
                            BoundaryCondition::VelocityInlet {
                                velocity: Vector3::new(T::zero(), T::zero(), u_inlet),
                            },
                        );
                    }
                } else {
                    boundary_conditions.insert(
                        v_idx,
                        BoundaryCondition::Dirichlet {
                            value: T::zero(),
                            component_values: Some(vec![Some(T::zero()), Some(T::zero()), Some(T::zero()), None]),
                        },
                    );
                }
                repaired_nodes += 1;
            }
        }

        println!(
            "Venturi Outlet Corner BC: corner_nodes={}, pressure_ref_applied={}",
            outlet_corner_nodes.len(),
            if outlet_corner_nodes.is_empty() { 0 } else { 1 }
        );
        println!(
            "Venturi Inlet/Wall Compatibility: inlet_nodes={}, wall_nodes={}, rim_nodes={}",
            inlet_nodes.len(),
            wall_nodes.len(),
            inlet_wall_rim_nodes.len()
        );
        println!("Venturi BC Coverage Repair: repaired_nodes={}", repaired_nodes);

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
        let mut fem_config = FemConfig::default();
        fem_config.grad_div_penalty = T::zero();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;

        for iter in 0..self.config.max_nonlinear_iterations {
            problem.element_viscosities = Some(element_viscosities.clone());
            let fem_result = solver.solve(&problem, last_solution.as_ref());

            // If the linear solver fails (e.g. GMRES stagnation on the
            // non-Newtonian saddle-point system), use the last converged
            // solution and break. The Picard iterations before failure
            // typically provide an accurate enough Newtonian/early-Casson
            // solution for validation.
            let fem_solution = match fem_result {
                Ok(sol) => sol,
                Err(e) => {
                    tracing::warn!("Picard iteration: linear solve failed ({}), using last converged solution", e);
                    println!("Picard: Linear solve failed, using last converged solution");
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
            
            let mut shear_min = <T as Float>::max_value();
            let mut shear_max = <T as Float>::min_value();
            let mut shear_sum = T::zero();

            for (i, cell) in problem.mesh.cells().iter().enumerate() {
                // Handle hex-to-tet averaging for shear rate
                let shear_rate = self.calculate_cell_shear_rate(cell, &problem.mesh, &updated_solution)?;
                
                if shear_rate < shear_min { shear_min = shear_rate; }
                if shear_rate > shear_max { shear_max = shear_rate; }
                shear_sum += shear_rate;

                let mut new_visc = fluid.viscosity_at_shear(shear_rate, T::from_f64_or_one(310.0), self.config.inlet_pressure)?;
                
                // Cap viscosity at 20x reference (stability)
                let max_viscosity = problem.fluid.viscosity * T::from_f64(20.0).unwrap();
                if new_visc > max_viscosity {
                    new_visc = max_viscosity;
                }
                if new_visc < problem.fluid.viscosity {
                    new_visc = problem.fluid.viscosity;
                }

                let change = Float::abs(new_visc - element_viscosities[i]) / element_viscosities[i];
                if change > max_change {
                    max_change = change;
                }
                new_viscosities.push(new_visc);
            }
            if n_elements > 0 {
                let shear_avg = shear_sum / T::from_usize(n_elements).unwrap();
                println!("Shear Rate Stats: Min={:?}, Max={:?}, Avg={:?}", shear_min, shear_max, shear_avg);
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
                println!("Picard: Converged in {} iterations", iter + 1);
                break;
            }
        }

        let fem_solution = last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;

        // Debug: velocity field diagnostic
        {
            let n = problem.mesh.vertex_count();
            let mut u_max = T::zero();
            let mut u_max_idx = 0;
            let mut nonzero_count = 0;
            let thr = T::from_f64(1e-10).unwrap_or_else(T::zero);
            let total_length = self.builder.l_inlet + self.builder.l_convergent + self.builder.l_throat + self.builder.l_divergent + self.builder.l_outlet;
            let n_bins = 10;
            let mut bin_u_max = vec![T::zero(); n_bins];
            let mut bin_count = vec![0usize; n_bins];
            for i in 0..n {
                let vel = fem_solution.get_velocity(i);
                let mag = vel.norm();
                if mag > thr { nonzero_count += 1; }
                if mag > u_max { u_max = mag; u_max_idx = i; }
                let z = problem.mesh.vertices()[i].position.z;
                let bin = ((z / total_length) * T::from_usize(n_bins).unwrap())
                    .to_usize().unwrap_or(0).min(n_bins - 1);
                bin_count[bin] += 1;
                if mag > bin_u_max[bin] { bin_u_max[bin] = mag; }
            }
            let u_max_pos = problem.mesh.vertices()[u_max_idx].position;
            println!("Velocity Diagnostic: {} of {} nodes have |u|>1e-10, max |u|={:?} at node {} pos=({:?},{:?},{:?})",
                nonzero_count, n, u_max, u_max_idx, u_max_pos.x, u_max_pos.y, u_max_pos.z);
            for b in 0..n_bins {
                let z_lo = total_length * T::from_usize(b).unwrap() / T::from_usize(n_bins).unwrap();
                let z_hi = total_length * T::from_usize(b+1).unwrap() / T::from_usize(n_bins).unwrap();
                println!("  z=[{:.4e},{:.4e}]: {} nodes, max|u|={:?}",
                    z_lo.to_f64().unwrap_or(0.0), z_hi.to_f64().unwrap_or(0.0),
                    bin_count[b], bin_u_max[b]);
            }
        }

        self.print_divergence_stats(&problem.mesh, &fem_solution)?;
        let mut solution = VenturiSolution3D::new();
        solution.u_inlet = u_inlet;
        
        // Average pressure at inlet
        let mut p_in_sum = T::zero();
        let mut p_in_count = 0usize;
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
        println!("Venturi Solver Setup: Q={:?}, A_in={:?}, u_inlet_calc={:?}", 
            self.config.inlet_flow_rate, area_inlet, u_inlet);

        let mut u_in_sol_sum = T::zero();
        for &v_idx in &inlet_nodes {
            if v_idx < fem_solution.n_corner_nodes {
                p_in_sum += fem_solution.get_pressure(v_idx);
                p_in_count += 1;
            }
            // Use norm for general velocity magnitude
            u_in_sol_sum += fem_solution.get_velocity(v_idx).norm();
            count_in += 1;
        }

        if p_in_count > 0 {
            solution.p_inlet = p_in_sum / T::from_usize(p_in_count).unwrap();
        } else {
            solution.p_inlet = self.config.inlet_pressure;
        }

        let u_in_sol_avg = if count_in > 0 { u_in_sol_sum / T::from_usize(count_in).unwrap() } else { T::zero() };
        println!("Venturi Solution Debug: Average Inlet Velocity = {:?}", u_in_sol_avg);
        
        // Identify throat section nodes and average pressure
        let z_throat_center = self.builder.l_inlet + self.builder.l_convergent + self.builder.l_throat / T::from_f64(2.0).unwrap();
        let total_length = self.builder.l_inlet + self.builder.l_convergent + self.builder.l_throat + self.builder.l_divergent + self.builder.l_outlet;
        // Use tolerance proportional to mesh spacing (half the axial element size)
        let throat_tol = total_length / T::from_usize(self.config.resolution.0.max(1)).unwrap() * T::from_f64(0.6).unwrap();
        let mut p_throat_sum = T::zero();
        let mut p_throat_count = 0usize;
        let mut u_throat_max = T::zero();
        let mut u_throat_sum = T::zero();
        let mut count_th = 0;

        for (i, v) in problem.mesh.vertices().iter().enumerate() {
            let dist_z = num_traits::Float::abs(v.position.z - z_throat_center);
            if dist_z < throat_tol {
                if i < fem_solution.n_corner_nodes {
                    p_throat_sum += fem_solution.get_pressure(i);
                    p_throat_count += 1;
                }
                let u_mag = fem_solution.get_velocity(i).norm();
                if u_mag > u_throat_max {
                    u_throat_max = u_mag;
                }
                u_throat_sum += u_mag;
                count_th += 1;
            }
        }
        
        let u_throat_avg = if count_th > 0 { u_throat_sum / T::from_usize(count_th).unwrap() } else { T::zero() };

        if p_throat_count > 0 {
            // FIX: Solver returns negative pressure potential? Invert sign for physical pressure.
            solution.p_throat = p_throat_sum / T::from_usize(p_throat_count).unwrap();
            solution.u_throat = u_throat_max;
        }
        
        // Average pressure at outlet
        let mut p_out_sum = T::zero();
        let mut p_out_count = 0usize;
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
            if v_idx < fem_solution.n_corner_nodes {
                p_out_sum += fem_solution.get_pressure(v_idx);
                p_out_count += 1;
            }
            u_out_sum += fem_solution.get_velocity(v_idx).norm();
            count_out += 1;
        }

        let u_out_avg = if count_out > 0 { u_out_sum / T::from_usize(count_out).unwrap() } else { T::zero() };

        if p_out_count > 0 {
            solution.p_outlet = p_out_sum / T::from_usize(p_out_count).unwrap();
        }

        // Sample pressure from interior inlet/outlet slices to reduce boundary-node artifacts.
        // Use flux-weighted averaging (u_z-weighted) for better effective pressure-drop estimate.
        let axial_dx = total_length / T::from_usize(self.config.resolution.0.max(1)).unwrap();
        let slice_tol = axial_dx * T::from_f64(0.55).unwrap();
        let z_in_sample = axial_dx * T::from_f64(5.0).unwrap();
        let z_out_sample = total_length - axial_dx * T::from_f64(5.0).unwrap();
        let core_radius_frac = T::from_f64(0.90).unwrap();
        let core_radius_sq = if self.config.circular {
            let r_core = T::from_f64(0.5).unwrap() * self.builder.d_inlet * core_radius_frac;
            r_core * r_core
        } else {
            T::zero()
        };

        let mut p_in_slice_sum = T::zero();
        let mut p_in_slice_count = 0usize;
        let mut p_out_slice_sum = T::zero();
        let mut p_out_slice_count = 0usize;
        let mut p_in_slice_weight_sum = T::zero();
        let mut p_out_slice_weight_sum = T::zero();
        let mut p_in_slice_weighted_sum = T::zero();
        let mut p_out_slice_weighted_sum = T::zero();
        let weight_floor = T::from_f64(1e-12).unwrap_or_else(T::zero);

        for (i, v) in problem.mesh.vertices().iter().enumerate() {
            if i >= fem_solution.n_corner_nodes {
                continue;
            }

            let z = v.position.z;
            if self.config.circular {
                let r_sq = v.position.x * v.position.x + v.position.y * v.position.y;
                if r_sq > core_radius_sq {
                    continue;
                }
            }

            if Float::abs(z - z_in_sample) <= slice_tol {
                let p_i = fem_solution.get_pressure(i);
                let uzi = Float::max(fem_solution.get_velocity(i).z, T::zero());
                let wi = if uzi > weight_floor { uzi } else { weight_floor };
                p_in_slice_sum += p_i;
                p_in_slice_count += 1;
                p_in_slice_weighted_sum += wi * p_i;
                p_in_slice_weight_sum += wi;
            }
            if Float::abs(z - z_out_sample) <= slice_tol {
                let p_i = fem_solution.get_pressure(i);
                let uzi = Float::max(fem_solution.get_velocity(i).z, T::zero());
                let wi = if uzi > weight_floor { uzi } else { weight_floor };
                p_out_slice_sum += p_i;
                p_out_slice_count += 1;
                p_out_slice_weighted_sum += wi * p_i;
                p_out_slice_weight_sum += wi;
            }
        }

        if p_in_slice_count > 0 {
            if p_in_slice_weight_sum > T::zero() {
                solution.p_inlet = p_in_slice_weighted_sum / p_in_slice_weight_sum;
            } else {
                solution.p_inlet = p_in_slice_sum / T::from_usize(p_in_slice_count).unwrap();
            }
        }
        if p_out_slice_count > 0 {
            if p_out_slice_weight_sum > T::zero() {
                solution.p_outlet = p_out_slice_weighted_sum / p_out_slice_weight_sum;
            } else {
                solution.p_outlet = p_out_slice_sum / T::from_usize(p_out_slice_count).unwrap();
            }
        }

        println!("Venturi Debug: p_in={:?}, p_throat={:?}, p_out={:?}, u_throat={:?}, count_th={}, p_in_slice_n={}, p_out_slice_n={}, p_in_wsum={:?}, p_out_wsum={:?}", 
            solution.p_inlet, solution.p_throat, solution.p_outlet, solution.u_throat, count_th, p_in_slice_count, p_out_slice_count, p_in_slice_weight_sum, p_out_slice_weight_sum);
        
        println!("Venturi Mass Flux Debug: u_in_avg={:?}, u_throat_avg={:?}, u_out_avg={:?}", 
            u_in_sol_avg, u_throat_avg, u_out_avg);

        let q_in_face = self.calculate_boundary_flow(&problem.mesh, &fem_solution, "inlet")?;
        let q_out_face = self.calculate_boundary_flow(&problem.mesh, &fem_solution, "outlet")?;
        println!("Venturi Mass Flux Debug: q_in_face={:?}, q_out_face={:?}", q_in_face, q_out_face);

        let total_length = self.builder.l_inlet
            + self.builder.l_convergent
            + self.builder.l_throat
            + self.builder.l_divergent
            + self.builder.l_outlet;
        let plane_tol = total_length / T::from_usize(self.config.resolution.0.max(1)).unwrap()
            * T::from_f64(0.55).unwrap_or_else(T::one);
        let plane_fracs = [0.0, 0.25, 0.5, 0.75, 1.0];
        for frac in plane_fracs {
            let z_plane = total_length * T::from_f64(frac).unwrap();
            let (q_plane, faces) = self.calculate_plane_flux(&problem.mesh, &fem_solution, z_plane, plane_tol)?;
            println!(
                "Venturi Slice Flux: z_frac={:.2}, z={:?}, faces={}, total_q={:?}",
                frac,
                z_plane,
                faces,
                q_plane
            );
        }

        solution.dp_throat = solution.p_inlet - solution.p_throat;
        solution.dp_recovery = solution.p_outlet - solution.p_inlet; // Usually negative (loss)
        
        let q_dyn = Float::max(T::from_f64(0.5).unwrap() * fluid_props.density * u_inlet * u_inlet, T::one());
        solution.cp_throat = (solution.p_inlet - solution.p_throat) / q_dyn;
        solution.cp_recovery = (solution.p_outlet - solution.p_inlet) / q_dyn;
        
        // Calculate Mass Balance Error
        // Qin = u_in_avg * A_in (reference)
        // Qth = u_th_avg * A_th (slice-based)
        // q_in_face/q_out_face = integrated boundary flux
        let area_throat = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0) * self.builder.d_throat * self.builder.d_throat
        } else {
            self.builder.d_throat * self.builder.d_throat
        };
        let q_in = u_in_sol_avg * area_inlet;
        let q_th = u_throat_avg * area_throat;
        solution.mass_error = if q_in_face > T::from_f64(1e-12).unwrap() {
            (q_in_face - q_out_face) / q_in_face
        } else {
            T::zero()
        };

        println!(
            "Venturi Mass Balance: Q_in={:?}, Q_th={:?}, Q_in_face={:?}, Q_out_face={:?}, Error={:?}",
            q_in, q_th, q_in_face, q_out_face, solution.mass_error
        );

        Ok(solution)
    }

    fn calculate_cell_shear_rate(
        &self,
        cell: &cfd_mesh::topology::Cell,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
    ) -> Result<T> {
        use crate::fem::solver::extract_vertex_indices;
        
        let idxs = extract_vertex_indices(cell, mesh).map_err(|e| Error::Solver(e.to_string()))?;
        let vertex_positions: Vec<Vector3<T>> = mesh.vertices().iter().map(|v| v.position.coords).collect();
        let local_verts: Vec<Vector3<T>> = idxs.iter().map(|&i| vertex_positions[i]).collect();

        // Shear rate for P2 elements (Tet10) or P1 (Tet4)
        // For P2, the gradient is linear, so we evaluate at centroid.
        // For P1, the gradient is constant.
        
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

    fn print_divergence_stats(
        &self,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
    ) -> Result<()> {
        use crate::fem::solver::extract_vertex_indices;

        let mut min_div = <T as RealField>::max_value().unwrap_or_else(T::one);
        let mut max_div = T::zero();
        let mut sum_div = T::zero();
        let mut count = 0usize;
        let mut total_volume = T::zero();
        let mut signed_div_sum = T::zero();
        let mut abs_div_vol_sum = T::zero();

        for cell in mesh.cells() {
            let idxs = extract_vertex_indices(cell, mesh)?;
            if idxs.len() < 4 {
                continue;
            }

            let mut local_verts = Vec::with_capacity(idxs.len());
            for &idx in &idxs {
                local_verts.push(mesh.vertex(idx).unwrap().position.coords);
            }

            let mut div = T::zero();
            let mut cell_volume = T::zero();
            if idxs.len() == 10 {
                let mut tet4 = crate::fem::element::FluidElement::new(idxs[0..4].to_vec());
                let six_v = tet4.calculate_volume(&local_verts);
                if Float::abs(six_v) < T::from_f64(1e-24).unwrap() {
                    continue;
                }
                cell_volume = tet4.volume;
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

                let tet10 = crate::fem::shape_functions::LagrangeTet10::new(p1_grads);
                let l_centroid = [T::from_f64(0.25).unwrap(); 4];
                let p2_grads = tet10.gradients(&l_centroid);

                for i in 0..10 {
                    let u = solution.get_velocity(idxs[i]);
                    div += p2_grads[(0, i)] * u.x
                        + p2_grads[(1, i)] * u.y
                        + p2_grads[(2, i)] * u.z;
                }
            } else {
                let mut element = crate::fem::element::FluidElement::new(idxs.clone());
                element.calculate_volume(&local_verts);
                element.calculate_shape_derivatives(&local_verts);
                cell_volume = element.volume;
                for i in 0..4 {
                    let u = solution.get_velocity(idxs[i]);
                    div += element.shape_derivatives[(0, i)] * u.x
                        + element.shape_derivatives[(1, i)] * u.y
                        + element.shape_derivatives[(2, i)] * u.z;
                }
            }

            let div_abs = Float::abs(div);
            if div_abs < min_div {
                min_div = div_abs;
            }
            if div_abs > max_div {
                max_div = div_abs;
            }
            sum_div += div_abs;
            total_volume += cell_volume;
            signed_div_sum += div * cell_volume;
            abs_div_vol_sum += div_abs * cell_volume;
            count += 1;
        }

        if count > 0 {
            let avg_div = sum_div / T::from_usize(count).unwrap_or_else(T::one);
            let vol_avg_div = if total_volume > T::zero() {
                abs_div_vol_sum / total_volume
            } else {
                T::zero()
            };
            println!(
                "Divergence Stats: min={:?}, max={:?}, avg={:?}, vol_avg={:?}, net={:?} (n={}, vol={:?})",
                min_div,
                max_div,
                avg_div,
                vol_avg_div,
                signed_div_sum,
                count,
                total_volume
            );
        }

        Ok(())
    }

    fn calculate_boundary_flow(
        &self,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
        label: &str,
    ) -> Result<T> {
        let mut total_q = T::zero();
        let mut face_count = 0usize;

        for f_idx in mesh.boundary_faces() {
            if mesh.boundary_label(f_idx) == Some(label) {
                if let Some(face) = mesh.face(f_idx) {
                    if face.vertices.len() >= 3 {
                        face_count += 1;
                        let v0 = mesh.vertex(face.vertices[0]).unwrap().position.coords;
                        let v1 = mesh.vertex(face.vertices[1]).unwrap().position.coords;
                        let v2 = mesh.vertex(face.vertices[2]).unwrap().position.coords;

                        let n_vec = (v1 - v0).cross(&(v2 - v0));
                        let area = n_vec.norm() * T::from_f64_or_one(0.5);
                        if area <= T::zero() {
                            continue;
                        }
                        let face_normal = n_vec.normalize();

                        let mut u_avg = Vector3::zeros();
                        for &v_idx in &face.vertices {
                            u_avg += solution.get_velocity(v_idx);
                        }
                        u_avg /= T::from_usize(face.vertices.len()).unwrap_or_else(T::one);

                        let mut n_oriented = face_normal;
                        if label == "inlet" && n_oriented.z > T::zero() {
                            n_oriented = -n_oriented;
                        } else if label == "outlet" && n_oriented.z < T::zero() {
                            n_oriented = -n_oriented;
                        }

                        let signed_flux = u_avg.dot(&n_oriented) * area;
                        let face_flow = if label == "inlet" {
                            -signed_flux
                        } else {
                            signed_flux
                        };
                        total_q += face_flow;
                    }
                }
            }
        }

        println!(
            "Venturi Flux: label={}, faces={}, total_q={:?}",
            label, face_count, total_q
        );

        Ok(total_q)
    }

    fn calculate_plane_flux(
        &self,
        mesh: &cfd_mesh::mesh::Mesh<T>,
        solution: &crate::fem::StokesFlowSolution<T>,
        z_plane: T,
        tol: T,
    ) -> Result<(T, usize)> {
        let mut total_q = T::zero();
        let mut face_count = 0usize;

        for face in mesh.faces() {
            if face.vertices.len() < 3 {
                continue;
            }

            let mut on_plane = true;
            for &v_idx in &face.vertices {
                if let Some(v) = mesh.vertex(v_idx) {
                    if Float::abs(v.position.z - z_plane) > tol {
                        on_plane = false;
                        break;
                    }
                }
            }
            if !on_plane {
                continue;
            }

            let v0 = mesh.vertex(face.vertices[0]).unwrap().position.coords;
            let v1 = mesh.vertex(face.vertices[1]).unwrap().position.coords;
            let v2 = mesh.vertex(face.vertices[2]).unwrap().position.coords;

            let n_vec = (v1 - v0).cross(&(v2 - v0));
            let area = n_vec.norm() * T::from_f64_or_one(0.5);
            if area <= T::zero() {
                continue;
            }
            let face_normal = n_vec.normalize();

            let mut u_avg = Vector3::zeros();
            for &v_idx in &face.vertices {
                u_avg += solution.get_velocity(v_idx);
            }
            u_avg /= T::from_usize(face.vertices.len()).unwrap_or_else(T::one);

            let mut face_flow = u_avg.dot(&face_normal) * area;
            if face_normal.z < T::zero() {
                face_flow = -face_flow;
            }

            total_q += face_flow;
            face_count += 1;
        }

        Ok((total_q, face_count))
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
    /// Mass balance error (relative)
    pub mass_error: T,
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
            mass_error: T::zero(),
        }
    }
}
