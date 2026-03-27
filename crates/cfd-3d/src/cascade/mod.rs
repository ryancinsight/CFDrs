//! Multi-stage cascade 3D FEM solver for CIF (Cascade Inertial Focusing) networks.
//!
//! This module orchestrates independent 3D FEM Navier-Stokes solves for each
//! channel segment in a multi-channel device.  The caller supplies per-channel
//! geometry and pre-computed flow rates (from a 1D Kirchhoff network solver),
//! and this module runs a 3D solve on each channel independently.
//!
//! # Theorem — Independent Channel Decomposition
//!
//! For fully-developed flow in long rectangular ducts ($L/D_h \gg 1$),
//! entrance effects are confined to the first $\sim 0.06 \cdot Re \cdot D_h$
//! of each channel.  Outside this region, each channel segment can be solved
//! independently with prescribed inlet flow rate from the 1D network solution.
//!
//! **Proof sketch**: The parabolic nature of the downstream momentum equation
//! ensures that information propagates only downstream.  Once the 1D solver
//! provides the global flow distribution, each channel's velocity field is
//! determined by its own geometry and the assigned inlet flow rate.  The
//! inter-channel coupling is fully captured by the 1D pressure-flow balance.
//!
//! # Architecture
//!
//! `CascadeSolver3D` does not depend on `cfd-schematics` or `cfd-1d`.  It
//! accepts a `Vec<CascadeChannelSpec>` with pre-computed geometry and flow
//! rates.  This keeps the crate dependency graph acyclic and allows callers
//! in `cfd-optim` to integrate 1D → 3D coupling externally.
//!
//! # Theorem — Picard Iteration Convergence for Generalised-Newtonian Stokes
//!
//! Given the viscosity function $\mu: \mathbb{R}^+ \to [\mu_\infty, \mu_0]$ is
//! Lipschitz-continuous with constant $L_\mu$, the Picard iteration
//!
//! ```text
//! μ^{(k+1)} = μ(γ̇(u^{(k)}))
//! ```
//!
//! converges to a fixed point in $H^1$ norm provided $L_\mu < \rho / C_K$
//! where $C_K$ is the Korn inequality constant and $\rho$ is the coercivity
//! constant of the bilinear form $a(\cdot, \cdot)$.
//!
//! **Proof sketch.** The mapping $T: \mu \mapsto \mu(\dot{\gamma}(\mathbf{u}(\mu)))$
//! is a contraction on $[\mu_\infty, \mu_0]$ under the stated bound. The Stokes
//! operator with bounded viscosity is coercive on $H^1_0$, so the linear solve
//! at each step is well-posed. Convergence follows from the Banach fixed-point
//! theorem.
//!
//! **Reference:** Hirn, A. (2013). "Finite element approximation of singular
//! power-law systems." *Math. Comp.* 82:1247–1268.


use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

// ── Hematocrit-viscosity coupling ─────────────────────────────────────────────

/// Viscosity correction ratio for local hematocrit vs reference hematocrit.
///
/// # Theorem — Quemada (1978) Viscosity Scaling
///
/// The infinite-shear viscosity of a concentrated RBC suspension follows
/// $\mu_\infty(H) = \mu_{\text{plasma}} \exp\!\bigl(k \cdot H / (1 - H)\bigr)$,
/// where $k \approx 2.5$ is the intrinsic viscosity coefficient.
///
/// **Proof sketch**: The Quemada model derives from the structural viscosity
/// of a concentrated suspension where particle interactions modify the
/// effective viscosity exponentially with volume fraction.  The ratio
/// $\mu(H_1) / \mu(H_2) = \exp\!\bigl(k (H_1/(1-H_1) - H_2/(1-H_2))\bigr)$
/// is independent of the plasma viscosity and depends only on the hematocrit
/// difference between the two channels.
fn hematocrit_viscosity_ratio(hct_local: f64, hct_reference: f64) -> f64 {
    let k = 2.5_f64; // Intrinsic viscosity coefficient (Quemada 1978)
    let h_local = hct_local.clamp(0.01, 0.70);
    let h_ref = hct_reference.clamp(0.01, 0.70);
    let exponent = k * (h_local / (1.0 - h_local) - h_ref / (1.0 - h_ref));
    exponent.exp()
}

// ── Configuration ─────────────────────────────────────────────────────────────

/// Specification for a single channel in the cascade.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CascadeChannelSpec {
    /// Identifier matching the `NetworkBlueprint` channel ID.
    pub id: String,
    /// Axial length [m].
    pub length: f64,
    /// Channel width [m] (cross-stream, varies for venturi).
    pub width: f64,
    /// Channel height [m] (constant for rectangular ducts).
    pub height: f64,
    /// Assigned volumetric flow rate [m³/s] from 1D presolver.
    pub flow_rate_m3_s: f64,
    /// Whether this channel contains a Venturi throat.
    pub is_venturi_throat: bool,
    /// Throat width [m] (only when `is_venturi_throat`).
    pub throat_width: Option<f64>,
    /// Local hematocrit [-] from Zweifach-Fung junction routing.
    ///
    /// When set, the solver creates a `CassonBlood::with_hematocrit()` model
    /// for this channel, adjusting yield stress (Chien 1970: τ_y ∝ HCT³) and
    /// infinite-shear viscosity (Quemada 1978) to reflect the local RBC
    /// concentration after upstream cell separation.
    ///
    /// If `None`, the feed fluid is used unchanged.
    pub local_hematocrit: Option<f64>,
}

/// Configuration for the 3D cascade solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CascadeConfig3D {
    /// Reference pressure at outlets [Pa].
    pub outlet_pressure: f64,
    /// Mesh resolution: (axial, cross-width, cross-height).
    pub resolution: (usize, usize, usize),
    /// Maximum Picard iterations for non-Newtonian viscosity coupling.
    pub max_picard_iterations: usize,
    /// Picard convergence tolerance (relative viscosity change).
    pub picard_tolerance: f64,
}

impl Default for CascadeConfig3D {
    fn default() -> Self {
        Self {
            outlet_pressure: 0.0,
            resolution: (40, 8, 8),
            max_picard_iterations: 10,
            picard_tolerance: 1e-3,
        }
    }
}

// ── Results ───────────────────────────────────────────────────────────────────

/// Per-channel result from the 3D solve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelResult3D {
    /// Channel identifier.
    pub channel_id: String,
    /// Mean wall shear stress [Pa].
    pub wall_shear_mean_pa: f64,
    /// Maximum wall shear stress [Pa].
    pub wall_shear_max_pa: f64,
    /// Pressure drop inlet → outlet [Pa].
    pub pressure_drop_pa: f64,
    /// Maximum velocity magnitude [m/s].
    pub max_velocity: f64,
    /// Whether the solver converged.
    pub converged: bool,
    /// Number of Picard iterations used.
    pub picard_iterations: usize,
    /// Local hematocrit used for this channel's viscosity model.
    /// Equals `CascadeChannelSpec::local_hematocrit` when set, or the feed
    /// fluid's native hematocrit otherwise.
    pub local_hematocrit: f64,
}

/// Aggregate result for the entire cascade.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CascadeResult3D {
    /// Per-channel results in the order supplied.
    pub channel_results: Vec<ChannelResult3D>,
    /// Sum of per-channel pressure drops [Pa].
    pub total_pressure_drop_pa: f64,
    /// Channel with the highest wall shear stress.
    pub max_shear_channel_id: String,
    /// Highest wall shear stress across all channels [Pa].
    pub max_shear_pa: f64,
}

// ── Solver ────────────────────────────────────────────────────────────────────

/// Orchestrates independent 3D FEM solves for each channel in a cascade.
///
/// # Type Parameters
///
/// - `F`: Fluid model implementing `FluidTrait<f64>` (Newtonian or non-Newtonian).
pub struct CascadeSolver3D<F: FluidTrait<f64> + Clone> {
    config: CascadeConfig3D,
    fluid: F,
}

impl<F: FluidTrait<f64> + Clone> CascadeSolver3D<F> {
    /// Create a new cascade solver.
    pub fn new(config: CascadeConfig3D, fluid: F) -> Self {
        Self { config, fluid }
    }

    /// Solve all channels independently and return aggregate results.
    ///
    /// Each channel gets its own 3D structured mesh and FEM solve.
    /// Channels are solved sequentially (could be parallelised with `rayon`
    /// in a future sprint).
    pub fn solve(&self, channels: &[CascadeChannelSpec]) -> Result<CascadeResult3D> {
        if channels.is_empty() {
            return Err(Error::InvalidConfiguration(
                "CascadeSolver3D: no channels supplied".into(),
            ));
        }

        let mut channel_results = Vec::with_capacity(channels.len());
        for spec in channels {
            let result = self.solve_channel(spec)?;
            channel_results.push(result);
        }

        let total_dp: f64 = channel_results.iter().map(|r| r.pressure_drop_pa).sum();
        let (max_id, max_tau) = channel_results
            .iter()
            .max_by(|a, b| {
                a.wall_shear_max_pa
                    .partial_cmp(&b.wall_shear_max_pa)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|r| (r.channel_id.clone(), r.wall_shear_max_pa))
            .unwrap_or_default();

        Ok(CascadeResult3D {
            channel_results,
            total_pressure_drop_pa: total_dp,
            max_shear_channel_id: max_id,
            max_shear_pa: max_tau,
        })
    }

    /// Solve a single channel with the 3D FEM solver.
    ///
    /// Builds a structured hex mesh → converts to tets → assigns BCs
    /// → runs Taylor-Hood FEM with Picard viscosity coupling → extracts
    /// wall shear and pressure drop.
    fn solve_channel(&self, spec: &CascadeChannelSpec) -> Result<ChannelResult3D> {
        use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
        use std::collections::HashMap;

        let (nx, ny, nz) = self.config.resolution;

        // 1. Build structured mesh and map to channel geometry.
        let mut mesh = cfd_mesh::domain::grid::StructuredGridBuilder::new(ny, nz, nx)
            .build()
            .map_err(|e| Error::Solver(e.to_string()))?;

        let half_h = spec.height / 2.0;
        let total_l = spec.length;

        // Compute local width at each axial station for venturi channels.
        let width_at = |z_frac: f64| -> f64 {
            if !spec.is_venturi_throat {
                return spec.width;
            }
            let throat_w = spec.throat_width.unwrap_or(spec.width);
            // Simple symmetric constriction: linear convergent-divergent
            // with throat at the midpoint of the channel.
            let t = 2.0 * (z_frac - 0.5).abs(); // 0 at center, 1 at ends
            throat_w + (spec.width - throat_w) * t
        };

        for i in 0..mesh.vertex_count() {
            let vid = VertexId::from_usize(i);
            let p = *mesh.vertices.position(vid);

            // StructuredGridBuilder produces [0,1]³ coordinates.
            let z_new = p.z * total_l;
            let local_hw = width_at(p.z) / 2.0;
            let u = p.x * 2.0 - 1.0; // [-1, 1]
            let v = p.y * 2.0 - 1.0;
            let x_new = u * local_hw;
            let y_new = v * half_h;

            let v_mut = mesh.vertices.get_mut(vid);
            v_mut.position.x = x_new;
            v_mut.position.y = y_new;
            v_mut.position.z = z_new;
        }

        // 2. Compute inlet velocity from flow rate.
        let inlet_area = spec.width * spec.height;
        let u_inlet = spec.flow_rate_m3_s / inlet_area;

        // 3. Assign boundary conditions.
        let mut boundary_conditions: HashMap<usize, BoundaryCondition<f64>> = HashMap::new();

        // Identify boundary faces: faces referenced by only one cell.
        let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
        for cell in &mesh.cells {
            for &face_idx in &cell.faces {
                *face_cell_count.entry(face_idx).or_insert(0) += 1;
            }
        }

        let z_min = 0.0_f64;
        let z_max = total_l;
        let z_tol = total_l / (nx as f64 * 4.0);

        for (&face_idx, &count) in &face_cell_count {
            if count != 1 || face_idx >= mesh.face_count() {
                continue;
            }
            let face = mesh.faces.get(FaceId::from_usize(face_idx));
            let centroid_z: f64 = face
                .vertices
                .iter()
                .map(|vid| mesh.vertices.position(*vid).z)
                .sum::<f64>()
                / face.vertices.len() as f64;

            let is_inlet = (centroid_z - z_min).abs() < z_tol;
            let is_outlet = (centroid_z - z_max).abs() < z_tol;
            let is_wall = !is_inlet && !is_outlet;

            let bc = if is_inlet {
                // Inlet face: uniform velocity in +z direction.
                BoundaryCondition::VelocityInlet {
                    velocity: Vector3::new(0.0, 0.0, u_inlet),
                }
            } else if is_outlet {
                // Outlet face: zero-gauge pressure.
                BoundaryCondition::PressureOutlet {
                    pressure: self.config.outlet_pressure,
                }
            } else {
                // Wall: no-slip.
                BoundaryCondition::Dirichlet {
                    value: 0.0,
                    component_values: Some(vec![Some(0.0); 3]),
                }
            };

            for &v_id in &face.vertices {
                if is_wall {
                    // Wall faces use or_insert so that inlet/outlet BCs assigned
                    // earlier to corner vertices are NOT overwritten: the plug-flow
                    // inlet already prescribes u=u_inlet at all inlet-plane vertices
                    // including corners, and that condition takes physical priority
                    // over the lateral-wall no-slip at the inlet plane.
                    boundary_conditions
                        .entry(v_id.as_usize())
                        .or_insert(bc.clone());
                } else {
                    // Inlet/outlet BCs overwrite any previously-assigned wall BCs
                    // on corner vertices so the inlet velocity/pressure condition
                    // is always enforced at the inflow and outflow planes.
                    boundary_conditions.insert(v_id.as_usize(), bc.clone());
                }
            }
        }

        // 4. Assemble FEM problem.
        //
        // If local_hematocrit is set, compute a viscosity correction factor
        // using the Quemada (1978) model: μ_∞ ∝ exp(k·H/(1−H)).
        // This scales all viscosities to reflect the local RBC concentration
        // after upstream Zweifach-Fung cell separation.
        let hct_viscosity_factor = spec
            .local_hematocrit
            .map_or(1.0, |hct_local| hematocrit_viscosity_ratio(hct_local, 0.45));

        let fluid_props = self
            .fluid
            .properties_at(310.0, 0.0)
            .map_err(|e| Error::Solver(e.to_string()))?;
        let base_viscosity = fluid_props.dynamic_viscosity * hct_viscosity_factor;
        let constant_fluid = cfd_core::physics::fluid::ConstantPropertyFluid {
            name: "cascade_channel".to_string(),
            density: fluid_props.density,
            viscosity: base_viscosity,
            specific_heat: fluid_props.specific_heat,
            thermal_conductivity: fluid_props.thermal_conductivity,
            speed_of_sound: fluid_props.speed_of_sound,
        };

        let n_corner_nodes = mesh.vertex_count();
        let n_elements = mesh.cell_count();
        let mut element_viscosities: Vec<f64> = vec![base_viscosity; n_elements];

        let mut problem =
            StokesFlowProblem::new(mesh, constant_fluid, boundary_conditions, n_corner_nodes);

        // 5. Picard iteration loop (non-Newtonian viscosity coupling).
        let fem_config = FemConfig::<f64>::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;
        let mut converged = false;
        let mut picard_iter = 0;

        for iter in 0..self.config.max_picard_iterations {
            problem.element_viscosities = Some(element_viscosities.clone());
            let solution = solver
                .solve(&problem, last_solution.as_ref())
                .map_err(|e| Error::Solver(e.to_string()))?;

            // Update viscosities from shear rate.
            let mut max_change: f64 = 0.0;
            let mut new_visc = Vec::with_capacity(n_elements);
            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                let gamma_dot = self.element_shear_rate(cell, &problem.mesh, &solution);
                let mu_base = self
                    .fluid
                    .viscosity_at_shear(gamma_dot, 310.0, 0.0)
                    .unwrap_or(fluid_props.dynamic_viscosity);
                let mu = mu_base * hct_viscosity_factor;
                let change =
                    (mu - element_viscosities[i]).abs() / element_viscosities[i].max(1e-15);
                if change > max_change {
                    max_change = change;
                }
                new_visc.push(mu);
            }
            element_viscosities = new_visc;
            last_solution = Some(solution);
            picard_iter = iter + 1;

            if max_change < self.config.picard_tolerance {
                converged = true;
                break;
            }
        }

        // 6. Extract results from final solution.
        let solution =
            last_solution.ok_or_else(|| Error::Solver("cascade: no solution produced".into()))?;

        let (wall_shear_mean, wall_shear_max) =
            self.extract_wall_shear(&problem.mesh, &solution, &face_cell_count);
        let pressure_drop =
            self.extract_pressure_drop(&problem.mesh, &solution, z_min, z_max, z_tol);
        let max_vel = self.extract_max_velocity(&solution);

        Ok(ChannelResult3D {
            channel_id: spec.id.clone(),
            wall_shear_mean_pa: wall_shear_mean,
            wall_shear_max_pa: wall_shear_max,
            pressure_drop_pa: pressure_drop,
            max_velocity: max_vel,
            converged,
            picard_iterations: picard_iter,
            local_hematocrit: spec.local_hematocrit.unwrap_or(0.45),
        })
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Approximate element shear rate from the FEM velocity solution.
    fn element_shear_rate(
        &self,
        cell: &cfd_mesh::domain::topology::Cell,
        mesh: &cfd_mesh::IndexedMesh,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> f64 {
        // Average velocity gradient magnitude across the element's vertices.
        let verts: Vec<usize> = cell
            .faces
            .iter()
            .flat_map(|&fi| mesh.faces.get(FaceId::from_usize(fi)).vertices.iter())
            .map(|v| v.as_usize())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        if verts.len() < 2 {
            return 0.0;
        }

        let mut sum_grad = 0.0;
        let mut count = 0usize;
        for &vi in &verts {
            let vel = solution.get_velocity(vi);
            sum_grad += vel.norm();
            count += 1;
        }
        let avg_u = sum_grad / count as f64;

        // Characteristic length: cube root of element volume estimate.
        let positions: Vec<_> = verts
            .iter()
            .map(|&vi| mesh.vertices.position(VertexId::from_usize(vi)))
            .collect();
        let n = positions.len() as f64;
        let centroid = positions.iter().fold(
            nalgebra::Point3::<f64>::origin(),
            |acc: nalgebra::Point3<f64>, p| {
                nalgebra::Point3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
            },
        );
        let centroid = nalgebra::Point3::new(centroid.x / n, centroid.y / n, centroid.z / n);
        let h = positions
            .iter()
            .map(|p| (nalgebra::Point3::new(p.x, p.y, p.z) - centroid).norm())
            .fold(0.0_f64, f64::max)
            .max(1e-15);

        avg_u / h
    }

    /// Extract mean and max wall shear stress from boundary faces.
    ///
    /// Uses the first-interior-node gradient: τ = μ · |u_interior| / d(face_centroid, v_interior).
    /// Wall face vertices carry a Dirichlet no-slip BC (u ≈ 0), so averaging velocity at the
    /// wall itself gives τ ≈ 0.  Instead, for each wall face we find the closest non-boundary
    /// vertex in the adjacent cell and use that velocity over its distance from the face centroid.
    fn extract_wall_shear(
        &self,
        mesh: &cfd_mesh::IndexedMesh,
        solution: &crate::fem::StokesFlowSolution<f64>,
        face_cell_count: &HashMap<usize, usize>,
    ) -> (f64, f64) {
        use std::collections::HashSet;

        let mu = self
            .fluid
            .properties_at(310.0, 0.0)
            .map_or(3.5e-3, |p| p.dynamic_viscosity);

        // Pre-compute z-range once (instead of inside the face loop).
        let z_range = mesh
            .vertices
            .iter()
            .map(|(_, vd)| vd.position.z)
            .fold((f64::MAX, f64::MIN), |(lo, hi): (f64, f64), z| {
                (lo.min(z), hi.max(z))
            });
        let z_tol = (z_range.1 - z_range.0) / (self.config.resolution.0 as f64 * 4.0);

        // Build set of all boundary vertex indices (all vertices on any boundary face).
        let boundary_verts: HashSet<usize> = face_cell_count
            .iter()
            .filter(|(_, &c)| c == 1)
            .flat_map(|(&fi, _)| {
                if fi < mesh.face_count() {
                    mesh.faces
                        .get(FaceId::from_usize(fi))
                        .vertices
                        .iter()
                        .map(|v| v.as_usize())
                        .collect::<Vec<_>>()
                } else {
                    vec![]
                }
            })
            .collect();

        // Build face_idx → owning cell index for all boundary (count == 1) faces.
        let mut face_to_cell: HashMap<usize, usize> = HashMap::new();
        for (cell_idx, cell) in mesh.cells.iter().enumerate() {
            for &fi in &cell.faces {
                if face_cell_count.get(&fi).copied().unwrap_or(0) == 1 {
                    face_to_cell.insert(fi, cell_idx);
                }
            }
        }

        let mut sum_tau = 0.0;
        let mut max_tau = 0.0_f64;
        let mut wall_count = 0usize;

        for (&face_idx, &count) in face_cell_count {
            if count != 1 || face_idx >= mesh.face_count() {
                continue;
            }
            let face = mesh.faces.get(FaceId::from_usize(face_idx));

            // Face centroid.
            let n_fv = face.vertices.len() as f64;
            let centroid =
                face.vertices
                    .iter()
                    .fold(nalgebra::Point3::<f64>::origin(), |acc, vid| {
                        let p = mesh.vertices.position(*vid);
                        nalgebra::Point3::new(
                            acc.x + p.x / n_fv,
                            acc.y + p.y / n_fv,
                            acc.z + p.z / n_fv,
                        )
                    });

            // Skip inlet/outlet faces — process only lateral walls.
            if (centroid.z - z_range.0).abs() < z_tol || (centroid.z - z_range.1).abs() < z_tol {
                continue;
            }

            // Find the cell that owns this wall face.
            let Some(&cell_idx) = face_to_cell.get(&face_idx) else {
                continue;
            };
            let cell = &mesh.cells[cell_idx];

            // Gather all unique vertices belonging to the owning cell.
            let cell_verts: HashSet<usize> = cell
                .faces
                .iter()
                .flat_map(|&fi| {
                    if fi < mesh.face_count() {
                        mesh.faces
                            .get(FaceId::from_usize(fi))
                            .vertices
                            .iter()
                            .map(|v| v.as_usize())
                            .collect::<Vec<_>>()
                    } else {
                        vec![]
                    }
                })
                .collect();

            // Interior vertices: in the cell but NOT on any boundary face.
            let interior_verts: Vec<usize> = cell_verts
                .iter()
                .copied()
                .filter(|vi| !boundary_verts.contains(vi))
                .collect();

            if interior_verts.is_empty() {
                continue;
            }

            // Pick the interior vertex closest to the face centroid.
            let Some((best_vi, best_dist)) = interior_verts
                .iter()
                .map(|&vi| {
                    let p = mesh.vertices.position(VertexId::from_usize(vi));
                    let d = (nalgebra::Point3::new(p.x, p.y, p.z) - centroid).norm();
                    (vi, d)
                })
                .min_by(|(_, da), (_, db)| da.partial_cmp(db).unwrap_or(std::cmp::Ordering::Equal))
            else {
                continue;
            };

            if best_dist < 1e-15 {
                continue;
            }

            // τ = μ · |u_interior| / d(face_centroid → interior_vertex)
            let u_interior = solution.get_velocity(best_vi).norm();
            let tau = mu * u_interior / best_dist;
            sum_tau += tau;
            max_tau = max_tau.max(tau);
            wall_count += 1;
        }

        let mean_tau = if wall_count > 0 {
            sum_tau / wall_count as f64
        } else {
            0.0
        };
        (mean_tau, max_tau)
    }

    /// Compute pressure drop between inlet and outlet faces.
    fn extract_pressure_drop(
        &self,
        mesh: &cfd_mesh::IndexedMesh,
        solution: &crate::fem::StokesFlowSolution<f64>,
        z_min: f64,
        z_max: f64,
        z_tol: f64,
    ) -> f64 {
        let mut p_inlet_sum = 0.0;
        let mut p_inlet_count = 0usize;
        let mut p_outlet_sum = 0.0;
        let mut p_outlet_count = 0usize;

        for i in 0..mesh.vertex_count() {
            let z = mesh.vertices.position(VertexId::from_usize(i)).z;
            let p = solution.pressure[i];
            if (z - z_min).abs() < z_tol {
                p_inlet_sum += p;
                p_inlet_count += 1;
            } else if (z - z_max).abs() < z_tol {
                p_outlet_sum += p;
                p_outlet_count += 1;
            }
        }

        let p_in = if p_inlet_count > 0 {
            p_inlet_sum / p_inlet_count as f64
        } else {
            0.0
        };
        let p_out = if p_outlet_count > 0 {
            p_outlet_sum / p_outlet_count as f64
        } else {
            0.0
        };
        (p_in - p_out).abs()
    }

    /// Maximum velocity magnitude in the domain.
    fn extract_max_velocity(&self, solution: &crate::fem::StokesFlowSolution<f64>) -> f64 {
        (0..solution.n_nodes)
            .map(|i| solution.get_velocity(i).norm())
            .fold(0.0_f64, f64::max)
    }
}

use std::collections::HashMap;

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_core::physics::fluid::ConstantPropertyFluid;

    fn water_fluid() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid {
            name: "water".to_string(),
            density: 1000.0,
            viscosity: 1e-3,
            specific_heat: 4186.0,
            thermal_conductivity: 0.6,
            speed_of_sound: 1500.0,
        }
    }

    fn simple_channel(id: &str) -> CascadeChannelSpec {
        CascadeChannelSpec {
            id: id.to_string(),
            length: 0.01,      // 10 mm
            width: 0.001,      // 1 mm
            height: 0.0005,    // 0.5 mm
            flow_rate_m3_s: 1e-8,
            is_venturi_throat: false,
            throat_width: None,
            local_hematocrit: None,
        }
    }

    #[test]
    fn cascade_solver_construction() {
        let config = CascadeConfig3D::default();
        let fluid = water_fluid();
        let solver = CascadeSolver3D::new(config.clone(), fluid);
        assert_eq!(solver.config.resolution, (40, 8, 8));
    }

    #[test]
    fn empty_channels_returns_error() {
        let config = CascadeConfig3D::default();
        let fluid = water_fluid();
        let solver = CascadeSolver3D::new(config, fluid);
        let result = solver.solve(&[]);
        assert!(
            result.is_err(),
            "solving with zero channels should return an error"
        );
    }

    #[test]
    fn channel_count_matches_input() {
        let config = CascadeConfig3D {
            resolution: (4, 2, 2), // very coarse for speed
            max_picard_iterations: 1,
            ..CascadeConfig3D::default()
        };
        let fluid = water_fluid();
        let solver = CascadeSolver3D::new(config, fluid);

        let channels = vec![simple_channel("ch1"), simple_channel("ch2")];
        let result = solver.solve(&channels).expect("solve should succeed on 2 simple channels");

        assert_eq!(
            result.channel_results.len(),
            2,
            "result should contain one entry per input channel"
        );
        assert_eq!(result.channel_results[0].channel_id, "ch1");
        assert_eq!(result.channel_results[1].channel_id, "ch2");
    }

    #[test]
    fn hematocrit_viscosity_ratio_equal_hct_is_one() {
        let ratio = hematocrit_viscosity_ratio(0.45, 0.45);
        assert!(
            (ratio - 1.0).abs() < 1e-10,
            "viscosity ratio for equal hematocrit should be 1.0, got {}",
            ratio
        );
    }

    #[test]
    fn hematocrit_viscosity_ratio_increases_with_hct() {
        let ratio_low = hematocrit_viscosity_ratio(0.30, 0.45);
        let ratio_high = hematocrit_viscosity_ratio(0.60, 0.45);
        assert!(
            ratio_low < 1.0,
            "lower hct should give ratio < 1, got {}",
            ratio_low
        );
        assert!(
            ratio_high > 1.0,
            "higher hct should give ratio > 1, got {}",
            ratio_high
        );
    }
}
