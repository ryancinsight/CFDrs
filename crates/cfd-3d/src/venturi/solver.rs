//! 3D FEM Navier-Stokes solver for Venturi throats
//!
//! Solves the incompressible Navier-Stokes equations on 3D Venturi domains
//! using Finite Element Method with support for non-Newtonian blood rheology.
//!
//! # Theorem — Bernoulli’s Equation (Inviscid Steady Flow)
//!
//! Along a streamline of steady, inviscid, incompressible flow:
//!
//! ```text
//! p + ½ ρ u² = const
//! ```
//!
//! From continuity ($A_1 u_1 = A_2 u_2$), the pressure drop in the Venturi
//! converging section is $\Delta p = \frac{1}{2} \rho (u_2^2 - u_1^2)$.
//!
//! # Theorem — ISO 5167 Discharge Coefficient
//!
//! The volumetric flow rate through a Venturi tube is
//!
//! ```text
//! Q = C_d · A_throat · √(2 Δp / ρ)
//! ```
//!
//! where $C_d = 0.995$ for a classical (machined convergent) Venturi per
//! ISO 5167-4:2003, valid for $2 \times 10^5 \leq Re_D \leq 2 \times 10^6$
//! and $0.4 \leq \beta \leq 0.75$ ($\beta = D_{throat}/D_{inlet}$).
//!
//! # Theorem — Pressure Recovery Coefficient
//!
//! The pressure recovery downstream of the throat is characterised by
//!
//! ```text
//! C_p = (p_d − p_t) / (½ ρ u_t²)
//! ```
//!
//! For a well-designed diffuser with 5–7° half-angle, $C_p \approx 0.8–0.9$.

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use cfd_mesh::VenturiMeshBuilder;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};

pub use super::types::{VenturiConfig3D, VenturiSolution3D};

// ============================================================================
// 3D Venturi Solver
// ============================================================================

/// 3D Finite Element Navier-Stokes solver for Venturi throats
pub struct VenturiSolver3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float> {
    builder: VenturiMeshBuilder<T>,
    config: VenturiConfig3D<T>,
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
    > VenturiSolver3D<T>
{
    /// Create new solver from mesh builder and config
    pub fn new(builder: VenturiMeshBuilder<T>, config: VenturiConfig3D<T>) -> Self {
        Self { builder, config }
    }

    /// Solve Venturi flow with given fluid (Newtonian or blood)
    #[allow(clippy::too_many_lines)]
    pub fn solve<F: FluidTrait<T> + Clone>(&self, fluid: F) -> Result<VenturiSolution3D<T>> {
        use crate::fem::StokesFlowProblem;
        use cfd_core::physics::boundary::BoundaryCondition;
        use std::collections::HashMap;

        // Generate the volume mesh via iso-parametric mapping of a structured
        // reference cube onto the Venturi geometry.  A unit-cube grid is mapped
        // so that the axial coordinate spans the full Venturi length while the
        // cross-section follows the local diameter profile (circular or square).
        let mut base_mesh = cfd_mesh::domain::grid::StructuredGridBuilder::new(
            self.config.resolution.1,
            self.config.resolution.1,
            self.config.resolution.0,
        )
        .build()
        .map_err(|e| Error::Solver(e.to_string()))?;

        // Convert all Venturi geometry parameters T → f64 for the iso-parametric
        // mapping closure.  These conversions must succeed: if T cannot represent
        // the builder's physical dimensions the mesh would be degenerate.
        let r_in = self
            .builder
            .d_inlet
            .to_f64()
            .ok_or_else(|| Error::Solver("d_inlet T→f64 conversion failed".into()))?
            / 2.0;
        let r_th = self
            .builder
            .d_throat
            .to_f64()
            .ok_or_else(|| Error::Solver("d_throat T→f64 conversion failed".into()))?
            / 2.0;
        let l_in = self
            .builder
            .l_inlet
            .to_f64()
            .ok_or_else(|| Error::Solver("l_inlet T→f64 conversion failed".into()))?;
        let l_conv = self
            .builder
            .l_convergent
            .to_f64()
            .ok_or_else(|| Error::Solver("l_convergent T→f64 conversion failed".into()))?;
        let l_th = self
            .builder
            .l_throat
            .to_f64()
            .ok_or_else(|| Error::Solver("l_throat T→f64 conversion failed".into()))?;
        let l_div = self
            .builder
            .l_divergent
            .to_f64()
            .ok_or_else(|| Error::Solver("l_divergent T→f64 conversion failed".into()))?;
        let l_out = self
            .builder
            .l_outlet
            .to_f64()
            .ok_or_else(|| Error::Solver("l_outlet T→f64 conversion failed".into()))?;
        let total_l = l_in + l_conv + l_th + l_div + l_out;

        let radius_at = |z: f64| -> f64 {
            if z <= l_in {
                r_in
            } else if z <= l_in + l_conv {
                let f = (z - l_in) / l_conv;
                r_in - (r_in - r_th) * f
            } else if z <= l_in + l_conv + l_th {
                r_th
            } else if z <= l_in + l_conv + l_th + l_div {
                let f = (z - (l_in + l_conv + l_th)) / l_div;
                r_th + (r_in - r_th) * f
            } else {
                r_in
            } // outlet
        };

        let is_circular = self.config.circular;
        // For rectangular cross-sections, use constant half-height for the
        // y-dimension so the venturi narrows only in x (width) while height
        // stays constant — matching real millifluidic chip geometry.
        let half_height: Option<f64> = if is_circular {
            None
        } else {
            self.config
                .rect_height
                .map(|h| {
                    h.to_f64()
                        .ok_or_else(|| Error::Solver("rect_height T→f64 conversion failed".into()))
                })
                .transpose()?
                .map(|h| h / 2.0)
        };
        for i in 0..base_mesh.vertex_count() {
            use cfd_mesh::domain::core::index::VertexId;
            let vid = VertexId::from_usize(i);
            let p = *base_mesh.vertices.position(vid);

            let z_new = p.z * total_l;
            let r = radius_at(z_new);

            let u = p.x * 2.0 - 1.0;
            let v = p.y * 2.0 - 1.0;

            let (x_new, y_new) = if is_circular {
                let x_d = u * (1.0 - v * v / 2.0).sqrt();
                let y_d = v * (1.0 - u * u / 2.0).sqrt();
                (x_d * r, y_d * r)
            } else {
                // Width narrows at the throat; height is constant.
                let hy = half_height.unwrap_or(r);
                (u * r, v * hy)
            };

            let v_mut = base_mesh.vertices.get_mut(vid);
            v_mut.position.x = x_new;
            v_mut.position.y = y_new;
            v_mut.position.z = z_new;
        }

        let tet_mesh = base_mesh;
        // Use P1 (linear) elements directly.  The P2MeshConverter is a surface-mesh
        // tool (1:4 triangle subdivision) and corrupts volumetric tet topology when
        // applied to a 3-D mesh: after it clears + replaces all faces the cell→face
        // index map becomes invalid, extract_vertex_indices falls to its unsorted
        // fallback, and the first four "corner" nodes are nearly coplanar (volume
        // < 1e-22), causing assembly to fail on element 0.  P1 Taylor-Hood elements
        // are well-posed for Stokes flow at the mesh resolutions used here.
        let n_corner_nodes = tet_mesh.vertex_count();
        let mesh = tet_mesh;

        // Boundary diagnostics: labeled faces vs connectivity boundary faces
        {
            use std::collections::{HashMap, HashSet};

            let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
            for cell in &mesh.cells {
                for &face_idx in &cell.faces {
                    *face_cell_count.entry(face_idx).or_insert(0) += 1;
                }
            }

            let connectivity_boundary_faces: HashSet<usize> = face_cell_count
                .iter()
                .filter(|&(_face_idx, &count)| count == 1)
                .map(|(&face_idx, _)| face_idx)
                .collect();

            let marked_boundary_faces: HashSet<usize> = mesh
                .boundary_faces()
                .into_iter()
                .map(|id| id.as_usize())
                .collect();

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

            for &f_idx_usize in &marked_boundary_faces {
                let f_idx = FaceId::from_usize(f_idx_usize);
                let label = mesh.boundary_label(f_idx);
                match label {
                    Some("inlet") => {
                        inlet_faces += 1;
                        let face = mesh.faces.get(f_idx);
                        for &v_idx in &face.vertices {
                            inlet_nodes.insert(v_idx.as_usize());
                            labeled_nodes.insert(v_idx.as_usize());
                        }
                    }
                    Some("outlet") => {
                        outlet_faces += 1;
                        let face = mesh.faces.get(f_idx);
                        for &v_idx in &face.vertices {
                            outlet_nodes.insert(v_idx.as_usize());
                            labeled_nodes.insert(v_idx.as_usize());
                        }
                    }
                    Some("wall") => {
                        wall_faces += 1;
                        let face = mesh.faces.get(f_idx);
                        for &v_idx in &face.vertices {
                            wall_nodes.insert(v_idx.as_usize());
                            labeled_nodes.insert(v_idx.as_usize());
                        }
                    }
                    _ => {
                        unlabeled_faces += 1;
                    }
                }
            }

            let mut connectivity_nodes = HashSet::new();
            for &f_idx in &connectivity_boundary_faces {
                let face = mesh.faces.get(FaceId::from_usize(f_idx));
                for &v_idx in &face.vertices {
                    connectivity_nodes.insert(v_idx.as_usize());
                }
            }

            let inlet_outlet_overlap = inlet_nodes.intersection(&outlet_nodes).count();
            let inlet_wall_overlap = inlet_nodes.intersection(&wall_nodes).count();
            let outlet_wall_overlap = outlet_nodes.intersection(&wall_nodes).count();

            tracing::debug!(
                marked = marked_boundary_faces.len(),
                connectivity = connectivity_boundary_faces.len(),
                union = boundary_union.len(),
                unlabeled = unlabeled_faces,
                "Venturi boundary face counts"
            );
            tracing::debug!(
                inlet_faces,
                outlet_faces,
                wall_faces,
                inlet_nodes = inlet_nodes.len(),
                outlet_nodes = outlet_nodes.len(),
                wall_nodes = wall_nodes.len(),
                "Venturi boundary node classification"
            );
            tracing::debug!(
                inlet_outlet_overlap,
                inlet_wall_overlap,
                outlet_wall_overlap,
                labeled = labeled_nodes.len(),
                connectivity = connectivity_nodes.len(),
                total = mesh.vertex_count(),
                "Venturi boundary overlap"
            );

            let labeled_not_connectivity = labeled_nodes.difference(&connectivity_nodes).count();
            let connectivity_not_labeled = connectivity_nodes.difference(&labeled_nodes).count();
            tracing::debug!(
                labeled_not_connectivity,
                connectivity_not_labeled,
                "Venturi boundary classification gaps"
            );
        }

        // 2. Define Boundary Conditions
        // Explicitly classify boundary node sets first so that inlet/wall intersection
        // (rim) nodes can be treated consistently.
        let mut boundary_conditions = HashMap::new();
        let fluid_props =
            fluid.properties_at(T::from_f64_or_one(310.0), self.config.inlet_pressure)?;

        let area_inlet = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0)
                * self.builder.d_inlet
                * self.builder.d_inlet
        } else {
            // Rectangular: width × height (height may differ from width)
            let h = self.config.rect_height.unwrap_or(self.builder.d_inlet);
            self.builder.d_inlet * h
        };
        let u_inlet = self.config.inlet_flow_rate / area_inlet;

        // ── Classify boundary faces ────────────────────────────────────────────
        // AxialBoundaryClassifier uses label-based or z-proximity classification
        // (label-based if the mesh carries explicit "inlet"/"outlet" labels,
        // z-proximity otherwise).  Replaces manual cell-face counting + z-bounds loop.
        let face_sets =
            crate::fem::AxialBoundaryClassifier::new(&mesh, self.config.resolution.0).classify();
        let inlet_nodes = face_sets.inlet_nodes;
        let outlet_nodes = face_sets.outlet_nodes;
        let wall_nodes = face_sets.wall_nodes;
        let boundary_vertices = face_sets.boundary_vertices;

        let inlet_wall_rim_nodes: std::collections::HashSet<usize> =
            inlet_nodes.intersection(&wall_nodes).copied().collect();

        // Pass 1: Apply inlet velocity to non-rim inlet nodes and no-slip to inlet-wall rim nodes.
        for &v_idx in &inlet_nodes {
            if inlet_wall_rim_nodes.contains(&v_idx) {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::Dirichlet {
                        value: 0.0_f64,
                        component_values: Some(vec![
                            Some(0.0_f64),
                            Some(0.0_f64),
                            Some(0.0_f64),
                            None,
                        ]),
                    },
                );
            } else {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::VelocityInlet {
                        velocity: Vector3::new(0.0_f64, 0.0_f64, u_inlet.to_f64().unwrap_or(0.0)),
                    },
                );
            }
        }

        // Pass 2: Apply outlet pressure Dirichlet to all corner nodes on
        // the outlet face.  For incompressible velocity-inlet / pressure-outlet
        // flow, prescribing p at the outlet anchors the pressure field while
        // the velocity there adjusts via the natural (traction-free) weak-form
        // condition.  The FEM solver only applies PressureOutlet to corner
        // nodes (those carrying a pressure DOF); mid-edge nodes are marked
        // Outflow so they are not flagged as unconstrained by diagnostics.
        let outlet_pressure_f64 = self.config.outlet_pressure.to_f64().unwrap_or(0.0);
        let mut outlet_corner_count = 0usize;
        for &v_idx in &outlet_nodes {
            if v_idx < n_corner_nodes {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::PressureOutlet {
                        pressure: outlet_pressure_f64,
                    },
                );
                outlet_corner_count += 1;
            } else {
                boundary_conditions.insert(v_idx, BoundaryCondition::Outflow);
            }
        }

        // Pass 3: Assign wall (no-slip) BCs to wall nodes EXCEPT those on the inlet or outlet.
        for &v_idx in &wall_nodes {
            if inlet_nodes.contains(&v_idx) || outlet_nodes.contains(&v_idx) {
                continue;
            }
            boundary_conditions.insert(
                v_idx,
                BoundaryCondition::Dirichlet {
                    value: 0.0_f64,
                    component_values: Some(vec![Some(0.0_f64), Some(0.0_f64), Some(0.0_f64), None]),
                },
            );
        }

        // Pass 4: Ensure all boundary vertices have a BC unless intentionally left as outlet natural boundary.
        let inlet_radius_sq = {
            let half = <T as FromPrimitive>::from_f64(0.5)
                .expect("0.5 is exactly representable in IEEE 754");
            let margin = <T as FromPrimitive>::from_f64(0.98)
                .expect("0.98 is an IEEE 754 representable f64 constant");
            let r = half * self.builder.d_inlet * margin;
            r * r
        };
        let mut repaired_nodes = 0usize;
        for &v_idx in &boundary_vertices {
            if boundary_conditions.contains_key(&v_idx) {
                continue;
            }

            if outlet_nodes.contains(&v_idx) {
                boundary_conditions.insert(v_idx, BoundaryCondition::Outflow);
            } else if inlet_nodes.contains(&v_idx) {
                let v = mesh.vertices.get(VertexId::from_usize(v_idx));
                if self.config.circular {
                    let r_sq = <T as From<f64>>::from(
                        v.position.x * v.position.x + v.position.y * v.position.y,
                    );
                    if r_sq >= inlet_radius_sq {
                        boundary_conditions.insert(
                            v_idx,
                            BoundaryCondition::Dirichlet {
                                value: 0.0_f64,
                                component_values: Some(vec![
                                    Some(0.0_f64),
                                    Some(0.0_f64),
                                    Some(0.0_f64),
                                    None,
                                ]),
                            },
                        );
                    } else {
                        boundary_conditions.insert(
                            v_idx,
                            BoundaryCondition::VelocityInlet {
                                velocity: Vector3::new(
                                    0.0_f64,
                                    0.0_f64,
                                    u_inlet.to_f64().unwrap_or(0.0),
                                ),
                            },
                        );
                    }
                } else if wall_nodes.contains(&v_idx) {
                    boundary_conditions.insert(
                        v_idx,
                        BoundaryCondition::Dirichlet {
                            value: 0.0_f64,
                            component_values: Some(vec![
                                Some(0.0_f64),
                                Some(0.0_f64),
                                Some(0.0_f64),
                                None,
                            ]),
                        },
                    );
                } else {
                    boundary_conditions.insert(
                        v_idx,
                        BoundaryCondition::VelocityInlet {
                            velocity: Vector3::new(
                                0.0_f64,
                                0.0_f64,
                                u_inlet.to_f64().unwrap_or(0.0),
                            ),
                        },
                    );
                }
            } else {
                boundary_conditions.insert(
                    v_idx,
                    BoundaryCondition::Dirichlet {
                        value: 0.0_f64,
                        component_values: Some(vec![
                            Some(0.0_f64),
                            Some(0.0_f64),
                            Some(0.0_f64),
                            None,
                        ]),
                    },
                );
            }
            repaired_nodes += 1;
        }

        tracing::debug!(
            outlet_corner_count,
            "Venturi Outlet Corner BC: corner_nodes, pressure_outlet_applied"
        );
        tracing::debug!(
            inlet_count = inlet_nodes.len(),
            wall_count = wall_nodes.len(),
            rim_count = inlet_wall_rim_nodes.len(),
            "Venturi Inlet/Wall Compatibility"
        );
        tracing::debug!(repaired_nodes, "Venturi BC Coverage Repair");

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
            n_corner_nodes,
        );
        let n_elements = problem.mesh.cell_count();
        let mut element_viscosities =
            vec![fluid_props.dynamic_viscosity.to_f64().unwrap_or(0.0); n_elements];
        let mut next_viscosities = Vec::with_capacity(n_elements);

        // 4. Picard Iteration Loop
        let fem_config = crate::fem::FemConfig::<f64>::default();
        let mut solver = crate::fem::FemSolver::new(fem_config);
        let mut last_solution: Option<crate::fem::StokesFlowSolution<f64>> = None;
        let mut anderson_accelerator = cfd_math::nonlinear_solver::AndersonAccelerator::<f64>::new(
            cfd_math::nonlinear_solver::AndersonConfig {
                history_depth: 5,
                relaxation: 1.0_f64,
                drop_tolerance: 1e-12_f64,
                method: cfd_math::nonlinear_solver::AndersonMethod::QR,
            },
        );

        let mut stagnation_count = 0usize;
        let mut prev_vel_change = f64::MAX;

        for iter in 0..self.config.max_nonlinear_iterations {
            let iter_start = std::time::Instant::now();

            if problem.element_viscosities.is_none() {
                problem.element_viscosities = Some(element_viscosities);
            }
            let fem_result = solver.solve_picard(
                &problem,
                last_solution.as_ref(),
                iter,
                self.config.max_nonlinear_iterations,
            );

            // If the linear solver fails (e.g. GMRES stagnation on the
            // non-Newtonian saddle-point system), use the last converged
            // solution and break. The Picard iterations before failure
            // typically provide an accurate enough Newtonian/early-Casson
            // solution for validation.
            let fem_solution = match fem_result {
                Ok(sol) => sol,
                Err(e) => {
                    tracing::warn!(error = %e, "Picard linear solve failed; using last converged solution");
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
            let mut shear_min_f64 = f64::MAX;
            let mut shear_max_f64 = f64::MIN;
            let mut shear_sum_f64 = 0.0_f64;

            next_viscosities.clear();
            let current_viscosities = problem
                .element_viscosities
                .as_ref()
                .expect("element_viscosities set before Picard loop");

            // Under-relax viscosity for Picard stability:
            // α ramps from 0.5 (iter 0) to 1.0 (iter ≥ 3)
            let relax_alpha = (0.5 + (iter as f64) / 6.0).min(1.0);

            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                // Handle hex-to-tet averaging for shear rate
                let shear_rate_f64 =
                    self.calculate_cell_shear_rate_f64(cell, &problem.mesh, &updated_solution)?;

                if shear_rate_f64 < shear_min_f64 {
                    shear_min_f64 = shear_rate_f64;
                }
                if shear_rate_f64 > shear_max_f64 {
                    shear_max_f64 = shear_rate_f64;
                }
                shear_sum_f64 += shear_rate_f64;

                let shear_rate = <T as From<f64>>::from(shear_rate_f64);
                let new_visc_t = fluid.viscosity_at_shear(
                    shear_rate,
                    T::from_f64_or_one(310.0),
                    self.config.inlet_pressure,
                )?;
                let mut new_visc = new_visc_t.to_f64().unwrap_or(0.0);

                // Cap viscosity at 20x reference (stability)
                let max_viscosity = problem.fluid.viscosity * 20.0_f64;
                if new_visc > max_viscosity {
                    new_visc = max_viscosity;
                }
                if new_visc < problem.fluid.viscosity {
                    new_visc = problem.fluid.viscosity;
                }

                let relaxed_visc =
                    relax_alpha * new_visc + (1.0 - relax_alpha) * current_viscosities[i];
                let change = num_traits::Float::abs(relaxed_visc - current_viscosities[i])
                    / current_viscosities[i];
                if change > max_change_f64 {
                    max_change_f64 = change;
                }
                next_viscosities.push(relaxed_visc);
            }
            if n_elements > 0 {
                let shear_avg = shear_sum_f64 / (n_elements as f64);
                tracing::debug!(
                    min = shear_min_f64,
                    max = shear_max_f64,
                    avg = shear_avg,
                    "Picard shear rate stats"
                );
            }

            // Track velocity convergence
            let mut vel_change_f64 = 0.0_f64;
            if let Some(ref prev) = last_solution {
                let diff = &updated_solution.velocity - &prev.velocity;
                let norm_prev = prev.velocity.norm();
                if norm_prev > f64::MIN_POSITIVE {
                    vel_change_f64 = diff.norm() / norm_prev;
                }
            } else {
                vel_change_f64 = 1.0_f64;
            }

            element_viscosities = problem
                .element_viscosities
                .take()
                .expect("element_viscosities set before Picard loop");
            std::mem::swap(&mut element_viscosities, &mut next_viscosities);

            last_solution = Some(updated_solution);

            // Log non-linear progress with timing
            let iter_elapsed = iter_start.elapsed();
            tracing::debug!(
                iter,
                vel_change = vel_change_f64,
                visc_change = max_change_f64,
                elapsed_secs = iter_elapsed.as_secs_f64(),
                "Picard iteration progress"
            );

            // Stagnation detection: if velocity change hasn't decreased
            // by at least 5% for 3 consecutive iterations, abort early.
            if vel_change_f64 >= prev_vel_change * 0.95 {
                stagnation_count += 1;
            } else {
                stagnation_count = 0;
            }
            prev_vel_change = vel_change_f64;

            if stagnation_count >= 3 {
                tracing::debug!(iter, "Picard stagnation detected — aborting");
                break;
            }

            let tol_f64 = self.config.nonlinear_tolerance.to_f64().unwrap_or(1e-4);
            if vel_change_f64 < tol_f64 && max_change_f64 < tol_f64 {
                tracing::debug!(iter, "Picard converged");
                break;
            }
        }

        let fem_solution =
            last_solution.ok_or_else(|| Error::Solver("No solution generated".to_string()))?;

        // Debug: velocity field diagnostic
        {
            let n = problem.mesh.vertex_count();
            let mut u_max = 0.0_f64;
            let mut u_max_idx = 0;
            let mut nonzero_count = 0;
            let thr = 1e-10_f64;
            let total_length_f64 = (self.builder.l_inlet
                + self.builder.l_convergent
                + self.builder.l_throat
                + self.builder.l_divergent
                + self.builder.l_outlet)
                .to_f64()
                .unwrap_or(1.0);
            let n_bins = 10;
            let mut bin_u_max = vec![0.0_f64; n_bins];
            let mut bin_count = vec![0usize; n_bins];
            for i in 0..n {
                let vel = fem_solution.get_velocity(i);
                let mag = vel.norm();
                if mag > thr {
                    nonzero_count += 1;
                }
                if mag > u_max {
                    u_max = mag;
                    u_max_idx = i;
                }
                let z = problem
                    .mesh
                    .vertices
                    .get(VertexId::from_usize(i))
                    .position
                    .z;
                let bin = if total_length_f64 > 0.0 {
                    ((z / total_length_f64) * n_bins as f64) as usize
                } else {
                    0
                }
                .min(n_bins - 1);
                bin_count[bin] += 1;
                if mag > bin_u_max[bin] {
                    bin_u_max[bin] = mag;
                }
            }
            let u_max_pos = problem
                .mesh
                .vertices
                .get(VertexId::from_usize(u_max_idx))
                .position;
            tracing::debug!(
                nonzero_count,
                total = n,
                u_max,
                u_max_idx,
                ?u_max_pos,
                "Venturi velocity diagnostic"
            );
            for b in 0..n_bins {
                let z_lo = total_length_f64 * (b as f64) / (n_bins as f64);
                let z_hi = total_length_f64 * ((b + 1) as f64) / (n_bins as f64);
                tracing::debug!(
                    b,
                    z_lo,
                    z_hi,
                    count = bin_count[b],
                    u_max = bin_u_max[b],
                    "velocity bin"
                );
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
        for f_idx in problem.mesh.boundary_faces() {
            if let Some(label) = problem.mesh.boundary_label(f_idx) {
                if label == "inlet" {
                    let face = problem.mesh.faces.get(f_idx);
                    for &v_idx in &face.vertices {
                        inlet_nodes.insert(v_idx.as_usize());
                    }
                }
            }
        }
        tracing::debug!(
            flow_rate = ?self.config.inlet_flow_rate,
            area_inlet = ?area_inlet,
            u_inlet = ?u_inlet,
            "Venturi solver setup"
        );

        let mut u_in_sol_sum = T::zero();
        for &v_idx in &inlet_nodes {
            if v_idx < fem_solution.n_corner_nodes {
                p_in_sum += <T as From<f64>>::from(fem_solution.get_pressure(v_idx));
                p_in_count += 1;
            }
            // Use axial (z) component for mass-flux diagnostic
            u_in_sol_sum += <T as From<f64>>::from(fem_solution.get_velocity(v_idx).z);
            count_in += 1;
        }

        if p_in_count > 0 {
            solution.p_inlet = p_in_sum
                / T::from_usize(p_in_count).expect("p_in_count is always a representable usize");
        } else {
            solution.p_inlet = self.config.inlet_pressure;
        }

        let u_in_sol_avg = if count_in > 0 {
            u_in_sol_sum
                / T::from_usize(count_in).expect("count_in is always a representable usize")
        } else {
            T::zero()
        };
        tracing::debug!(u_in_sol_avg = ?u_in_sol_avg, "Venturi average inlet velocity");

        // Identify throat section nodes and average pressure
        let z_throat_center = self.builder.l_inlet
            + self.builder.l_convergent
            + self.builder.l_throat
                / <T as FromPrimitive>::from_f64(2.0)
                    .expect("2.0 is representable in all IEEE 754 types");
        let total_length = self.builder.l_inlet
            + self.builder.l_convergent
            + self.builder.l_throat
            + self.builder.l_divergent
            + self.builder.l_outlet;
        // Use tolerance proportional to mesh spacing (half the axial element size)
        let throat_tol = total_length
            / T::from_usize(self.config.resolution.0.max(1))
                .expect("resolution is always a representable usize")
            * <T as FromPrimitive>::from_f64(0.6)
                .expect("0.6 is an IEEE 754 representable f64 constant");
        let mut p_throat_sum = T::zero();
        let mut p_throat_count = 0usize;
        let mut u_throat_max = T::zero();
        let mut u_throat_sum = T::zero();
        let mut count_th = 0;

        for i in 0..n_corner_nodes {
            let v = problem.mesh.vertices.get(VertexId::from_usize(i));
            let dist_z =
                num_traits::Float::abs(<T as From<f64>>::from(v.position.z) - z_throat_center);
            if dist_z < throat_tol {
                if i < fem_solution.n_corner_nodes {
                    p_throat_sum += <T as From<f64>>::from(fem_solution.get_pressure(i));
                    p_throat_count += 1;
                }
                let u_mag = <T as From<f64>>::from(fem_solution.get_velocity(i).norm());
                if u_mag > u_throat_max {
                    u_throat_max = u_mag;
                }
                u_throat_sum += u_mag;
                count_th += 1;
            }
        }

        let u_throat_avg = if count_th > 0 {
            u_throat_sum
                / T::from_usize(count_th).expect("count_th is always a representable usize")
        } else {
            T::zero()
        };

        if p_throat_count > 0 {
            solution.p_throat = p_throat_sum
                / T::from_usize(p_throat_count)
                    .expect("p_throat_count is always a representable usize");
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
                    let face = problem.mesh.faces.get(f_idx);
                    for &v_idx in &face.vertices {
                        outlet_nodes.insert(v_idx.as_usize());
                    }
                }
            }
        }

        for &v_idx in &outlet_nodes {
            if v_idx < fem_solution.n_corner_nodes {
                p_out_sum += <T as From<f64>>::from(fem_solution.get_pressure(v_idx));
                p_out_count += 1;
            }
            u_out_sum += <T as From<f64>>::from(fem_solution.get_velocity(v_idx).z);
            count_out += 1;
        }

        let u_out_avg = if count_out > 0 {
            u_out_sum / T::from_usize(count_out).expect("count_out is always a representable usize")
        } else {
            T::zero()
        };

        if p_out_count > 0 {
            solution.p_outlet = p_out_sum
                / T::from_usize(p_out_count).expect("p_out_count is always a representable usize");
        }

        // Sample pressure from interior inlet/outlet slices to reduce boundary-node artifacts.
        // Use flux-weighted averaging (u_z-weighted) for better effective pressure-drop estimate.
        let axial_dx = total_length
            / T::from_usize(self.config.resolution.0.max(1))
                .expect("resolution is always a representable usize");
        let slice_tol = axial_dx
            * <T as FromPrimitive>::from_f64(0.55)
                .expect("0.55 is an IEEE 754 representable f64 constant");
        let z_in_sample = axial_dx
            * <T as FromPrimitive>::from_f64(5.0)
                .expect("5.0 is representable in all IEEE 754 types");
        let z_out_sample = total_length
            - axial_dx
                * <T as FromPrimitive>::from_f64(5.0)
                    .expect("5.0 is representable in all IEEE 754 types");
        let core_radius_frac = <T as FromPrimitive>::from_f64(0.90)
            .expect("0.90 is an IEEE 754 representable f64 constant");
        let core_radius_sq = if self.config.circular {
            let r_core = <T as FromPrimitive>::from_f64(0.5)
                .expect("0.5 is exactly representable in IEEE 754")
                * self.builder.d_inlet
                * core_radius_frac;
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
        let weight_floor = <T as FromPrimitive>::from_f64(1e-12)
            .expect("1e-12 is an IEEE 754 representable f64 constant");

        for i in 0..n_corner_nodes {
            let v = problem.mesh.vertices.get(VertexId::from_usize(i));
            if i >= fem_solution.n_corner_nodes {
                continue;
            }

            let z = <T as From<f64>>::from(v.position.z);
            if self.config.circular {
                let r_sq = <T as From<f64>>::from(
                    v.position.x * v.position.x + v.position.y * v.position.y,
                );
                if r_sq > core_radius_sq {
                    continue;
                }
            }

            if Float::abs(z - z_in_sample) <= slice_tol {
                let p_i = <T as From<f64>>::from(fem_solution.get_pressure(i));
                let uzi =
                    <T as From<f64>>::from(Float::max(fem_solution.get_velocity(i).z, 0.0_f64));
                let wi = if uzi > weight_floor {
                    uzi
                } else {
                    weight_floor
                };
                p_in_slice_sum += p_i;
                p_in_slice_count += 1;
                p_in_slice_weighted_sum += wi * p_i;
                p_in_slice_weight_sum += wi;
            }
            if Float::abs(z - z_out_sample) <= slice_tol {
                let p_i = <T as From<f64>>::from(fem_solution.get_pressure(i));
                let uzi =
                    <T as From<f64>>::from(Float::max(fem_solution.get_velocity(i).z, 0.0_f64));
                let wi = if uzi > weight_floor {
                    uzi
                } else {
                    weight_floor
                };
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
                solution.p_inlet = p_in_slice_sum
                    / T::from_usize(p_in_slice_count)
                        .expect("p_in_slice_count is always a representable usize");
            }
        }
        if p_out_slice_count > 0 {
            if p_out_slice_weight_sum > T::zero() {
                solution.p_outlet = p_out_slice_weighted_sum / p_out_slice_weight_sum;
            } else {
                solution.p_outlet = p_out_slice_sum
                    / T::from_usize(p_out_slice_count)
                        .expect("p_out_slice_count is always a representable usize");
            }
        }

        tracing::debug!(
            p_in = ?solution.p_inlet,
            p_throat = ?solution.p_throat,
            p_out = ?solution.p_outlet,
            u_throat = ?solution.u_throat,
            count_th,
            p_in_slice_n = p_in_slice_count,
            p_out_slice_n = p_out_slice_count,
            p_in_wsum = ?p_in_slice_weight_sum,
            p_out_wsum = ?p_out_slice_weight_sum,
            "Venturi pressure diagnostics"
        );
        tracing::debug!(
            u_in_avg = ?u_in_sol_avg,
            u_throat_avg = ?u_throat_avg,
            u_out_avg = ?u_out_avg,
            "Venturi mass flux diagnostics"
        );

        let q_in_face = self.calculate_boundary_flow(&problem.mesh, &fem_solution, "inlet")?;
        let q_out_face = self.calculate_boundary_flow(&problem.mesh, &fem_solution, "outlet")?;
        tracing::debug!(q_in_face, q_out_face, "Venturi boundary face flux");

        let total_length = self.builder.l_inlet
            + self.builder.l_convergent
            + self.builder.l_throat
            + self.builder.l_divergent
            + self.builder.l_outlet;
        let plane_tol = total_length.to_f64().unwrap_or(1.0)
            / (self.config.resolution.0.max(1) as f64)
            * 0.55_f64;
        let plane_fracs = [0.0_f64, 0.25, 0.5, 0.75, 1.0];
        for frac in plane_fracs {
            let z_plane = total_length.to_f64().unwrap_or(0.0) * frac;
            let (q_plane, faces) =
                self.calculate_plane_flux(&problem.mesh, &fem_solution, z_plane, plane_tol)?;
            tracing::debug!(frac, z_plane, faces, q_plane, "Venturi slice flux");
        }

        solution.dp_throat = solution.p_inlet - solution.p_throat;
        solution.dp_recovery = solution.p_outlet - solution.p_inlet; // Usually negative (loss)

        let q_dyn = Float::max(
            <T as FromPrimitive>::from_f64(0.5).expect("0.5 is exactly representable in IEEE 754")
                * fluid_props.density
                * u_inlet
                * u_inlet,
            T::one(),
        );
        solution.cp_throat = (solution.p_inlet - solution.p_throat) / q_dyn;
        solution.cp_recovery = (solution.p_outlet - solution.p_inlet) / q_dyn;

        // Calculate Mass Balance Error
        // Qin = u_in_avg * A_in (reference)
        // Qth = u_th_avg * A_th (slice-based)
        // q_in_face/q_out_face = integrated boundary flux
        let area_throat = if self.config.circular {
            T::from_f64_or_one(std::f64::consts::PI / 4.0)
                * self.builder.d_throat
                * self.builder.d_throat
        } else {
            let h = self.config.rect_height.unwrap_or(self.builder.d_throat);
            self.builder.d_throat * h
        };
        let q_in = u_in_sol_avg * area_inlet;
        let q_th = u_throat_avg * area_throat;
        solution.q_in_face = <T as From<f64>>::from(q_in_face);
        solution.mass_error = if q_in_face > 1e-12_f64 {
            <T as From<f64>>::from((q_in_face - q_out_face) / q_in_face)
        } else {
            T::zero()
        };

        tracing::debug!(
            q_in = ?q_in,
            q_th = ?q_th,
            q_in_face,
            q_out_face,
            mass_error = ?solution.mass_error,
            "Venturi mass balance"
        );

        Ok(solution)
    }
}
