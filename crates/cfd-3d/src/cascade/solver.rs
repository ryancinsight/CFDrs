#![allow(clippy::uninlined_format_args)]
use std::collections::HashMap;

use aequitas::systems::si::quantities::{Area, Pressure, Velocity};
use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use leto::geometry::Vector3;

use super::types::{CascadeChannelSpec, CascadeConfig3D, CascadeResult3D, ChannelResult3D};

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
                     // Lower-bound at zero so plasma-dominant lanes are not artificially floored.
    let h_local = hct_local.clamp(0.0, 0.70);
    let h_ref = hct_reference.clamp(0.0, 0.70);
    let exponent = k * (h_local / (1.0 - h_local) - h_ref / (1.0 - h_ref));
    exponent.exp()
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
    /// Channels are solved sequentially; future parallel channel evaluation
    /// must route through the workspace Moirai execution provider.
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

        let total_dp = Pressure::from_base(
            channel_results
                .iter()
                .map(|r| r.pressure_drop_pa.into_base())
                .sum(),
        );
        let (max_id, max_tau) = channel_results
            .iter()
            .max_by(|a, b| {
                a.wall_shear_max_pa
                    .into_base()
                    .partial_cmp(&b.wall_shear_max_pa.into_base())
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

        let half_h = spec.height.into_base() / 2.0;
        let total_l = spec.length.into_base();

        // Compute local width at each axial station for venturi channels.
        let width_at = |z_frac: f64| -> f64 {
            if !spec.is_venturi_throat {
                return spec.width.into_base();
            }
            let throat_w = spec.throat_width.unwrap_or(spec.width).into_base();
            // Simple symmetric constriction: linear convergent-divergent
            // with throat at the midpoint of the channel.
            let t = 2.0 * (z_frac - 0.5).abs(); // 0 at center, 1 at ends
            throat_w + (spec.width.into_base() - throat_w) * t
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
        let inlet_area = Area::from_base(spec.width.into_base() * spec.height.into_base());
        let inlet_velocity: Velocity = spec.flow_rate_m3_s / inlet_area;
        let u_inlet = inlet_velocity.into_base();

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
                    pressure: self.config.outlet_pressure.into_base(),
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
        let base_viscosity = fluid_props.dynamic_viscosity.into_base() * hct_viscosity_factor;
        let constant_fluid = cfd_core::physics::fluid::ConstantPropertyFluid::new(
            "cascade_channel".to_string(),
            aequitas::systems::si::quantities::MassDensity::from_base(
                fluid_props.density.into_base(),
            ),
            aequitas::systems::si::quantities::DynamicViscosity::from_base(base_viscosity),
            aequitas::systems::si::quantities::SpecificHeatCapacity::from_base(
                fluid_props.specific_heat.into_base(),
            ),
            aequitas::systems::si::quantities::ThermalConductivity::from_base(
                fluid_props.thermal_conductivity.into_base(),
            ),
            aequitas::systems::si::quantities::Velocity::from_base(
                fluid_props.speed_of_sound.into_base(),
            ),
        );

        let n_corner_nodes = mesh.vertex_count();
        let n_elements = mesh.cell_count();
        let mut element_viscosities: Vec<f64> = vec![base_viscosity; n_elements];
        let mut next_element_viscosities: Vec<f64> = Vec::with_capacity(n_elements);

        let mut problem =
            StokesFlowProblem::new(mesh, constant_fluid, boundary_conditions, n_corner_nodes);

        // 5. Picard iteration loop (non-Newtonian viscosity coupling).
        let fem_config = FemConfig::<f64>::default();
        let mut solver = FemSolver::new(fem_config);
        let mut last_solution = None;
        let mut converged = false;
        let mut picard_iter = 0;

        for iter in 0..self.config.max_picard_iterations {
            problem.element_viscosities = Some(std::mem::replace(
                &mut element_viscosities,
                Vec::with_capacity(n_elements),
            ));
            let solution = solver
                .solve(&problem, last_solution.as_ref())
                .map_err(|e| Error::Solver(e.to_string()))?;

            let old_viscosities = problem
                .element_viscosities
                .take()
                .expect("element viscosities are assigned before Picard solve");

            // Update viscosities from shear rate.
            let mut max_change: f64 = 0.0;
            next_element_viscosities.clear();
            for (i, cell) in problem.mesh.cells.iter().enumerate() {
                let gamma_dot = self.element_shear_rate(cell, &problem.mesh, &solution);
                let mu_base = self
                    .fluid
                    .viscosity_at_shear(gamma_dot, 310.0, 0.0)
                    .map_or(fluid_props.dynamic_viscosity.into_base(), |value| {
                        value.into_base()
                    });
                let mu = mu_base * hct_viscosity_factor;
                let change = (mu - old_viscosities[i]).abs() / old_viscosities[i].max(1e-15);
                if change > max_change {
                    max_change = change;
                }
                next_element_viscosities.push(mu);
            }
            std::mem::swap(&mut element_viscosities, &mut next_element_viscosities);
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
            wall_shear_mean_pa: Pressure::from_base(wall_shear_mean),
            wall_shear_max_pa: Pressure::from_base(wall_shear_max),
            pressure_drop_pa: Pressure::from_base(pressure_drop),
            max_velocity: Velocity::from_base(max_vel),
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
        use std::collections::HashSet;

        // Average velocity gradient magnitude across the element's unique vertices.
        let mut unique_verts = HashSet::with_capacity(cell.faces.len() * 4);
        for &fi in &cell.faces {
            if fi >= mesh.face_count() {
                continue;
            }
            unique_verts.extend(
                mesh.faces
                    .get(FaceId::from_usize(fi))
                    .vertices
                    .iter()
                    .map(|v| v.as_usize()),
            );
        }
        if unique_verts.len() < 2 {
            return 0.0;
        }

        let mut sum_grad = 0.0;
        let mut positions = Vec::with_capacity(unique_verts.len());
        for &vi in &unique_verts {
            let vel = solution.get_velocity(vi);
            sum_grad += vel.norm();
            positions.push(*mesh.vertices.position(VertexId::from_usize(vi)));
        }
        let avg_u = sum_grad / unique_verts.len() as f64;

        // Characteristic length: cube root of element volume estimate.
        let n = positions.len() as f64;
        let centroid = positions.iter().fold(
            leto::Point3::<f64>::origin(),
            |acc: leto::Point3<f64>, p| leto::Point3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z),
        );
        let centroid = leto::Point3::new(centroid.x / n, centroid.y / n, centroid.z / n);
        let h = positions
            .iter()
            .map(|p| (leto::Point3::new(p.x, p.y, p.z) - centroid).norm())
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
            .map_or(3.5e-3, |p| p.dynamic_viscosity.into_base());

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
            let centroid = face
                .vertices
                .iter()
                .fold(leto::Point3::<f64>::origin(), |acc, vid| {
                    let p = mesh.vertices.position(*vid);
                    leto::Point3::new(acc.x + p.x / n_fv, acc.y + p.y / n_fv, acc.z + p.z / n_fv)
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
                    let d = (leto::Point3::new(p.x, p.y, p.z) - centroid).norm();
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

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::quantities::{Length, VolumetricFlowRate};
    use cfd_core::physics::fluid::ConstantPropertyFluid;

    fn water_fluid() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "water".to_string(),
            aequitas::systems::si::quantities::MassDensity::from_base(1000.0),
            aequitas::systems::si::quantities::DynamicViscosity::from_base(1e-3),
            aequitas::systems::si::quantities::SpecificHeatCapacity::from_base(4186.0),
            aequitas::systems::si::quantities::ThermalConductivity::from_base(0.6),
            aequitas::systems::si::quantities::Velocity::from_base(1500.0),
        )
    }

    fn simple_channel(id: &str) -> CascadeChannelSpec {
        CascadeChannelSpec {
            id: id.to_string(),
            length: Length::from_base(0.01),   // 10 mm
            width: Length::from_base(0.001),   // 1 mm
            height: Length::from_base(0.0005), // 0.5 mm
            flow_rate_m3_s: VolumetricFlowRate::from_base(1e-8),
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
            outlet_pressure: Pressure::from_base(0.0),
            resolution: (4, 2, 2), // very coarse for speed
            max_picard_iterations: 1,
            ..CascadeConfig3D::default()
        };
        let fluid = water_fluid();
        let solver = CascadeSolver3D::new(config, fluid);

        let channels = vec![simple_channel("ch1"), simple_channel("ch2")];
        let result = solver
            .solve(&channels)
            .expect("solve should succeed on 2 simple channels");

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

    #[test]
    fn hematocrit_viscosity_ratio_zero_local_hematocrit_matches_analytic_ratio() {
        let hct_reference = 0.45_f64;
        let k = 2.5_f64;
        let expected = (-k * hct_reference / (1.0 - hct_reference)).exp();
        let ratio = hematocrit_viscosity_ratio(0.0, hct_reference);

        assert!(
            (ratio - expected).abs() < 1e-15 * expected.max(1.0),
            "expected {expected}, got {ratio}"
        );
    }
}
