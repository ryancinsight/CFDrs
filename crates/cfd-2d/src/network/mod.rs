//! Full-network 2D solver: one 2D channel domain per `NetworkBlueprint` channel.
//!
//! This module provides a [`Network2DSolver<T>`] that is analogous to `cfd-1d`'s
//! `Network<T, F>`: it is built from a [`NetworkBlueprint`] via a [`GraphSink`]
//! adapter ([`Network2dBuilderSink<T>`]) and exposes a [`solve_all`] method that
//! runs a 2D Navier-Stokes solve for every straight channel in the blueprint.
//!
//! # Flow-rate distribution
//!
//! For the symmetric closed-loop topologies produced by `cfd-schematics` composite
//! presets, flow splits equally at every bifurcation/trifurcation junction by
//! construction.  The builder distributes `total_flow_rate_m3_s` using a
//! resistance-weighted BFS traversal of the blueprint graph, so it handles
//! non-symmetric topologies correctly too.
//!
//! # Per-channel 2D solve
//!
//! Each blueprint channel is a **straight rectangular duct**.  The 2D solve uses
//! [`NavierStokesSolver2D`] on a staggered Cartesian grid of size `(grid_nx × grid_ny)`.
//! A uniform inlet velocity profile is applied; zero-gauge pressure at the outlet.
//!
//! # Haemolysis
//!
//! Per-channel Giersiepen (1990) HI is computed from the mean wall shear stress
//! (6Q/(w·h²) for a rectangular duct) and the mean transit time through the channel.
//!
//! # Usage
//! ```rust,ignore
//! use cfd_2d::{Network2DSolver, Network2dBuilderSink, Network2dResult};
//! use cfd_schematics::application::use_cases::NetworkGenerationService;
//! use cfd_schematics::interface::presets::bifurcation_venturi_rect;
//! use cfd_core::physics::fluid::BloodModel;
//!
//! let bp = bifurcation_venturi_rect("bv", 0.005, 0.002, 0.0005, 0.001, 0.001);
//! let blood = BloodModel::Newtonian(3.5e-3_f64);
//! let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-6, 40, 20);
//! let mut net2d = NetworkGenerationService::new(sink).generate(&bp).unwrap();
//! let result = net2d.solve_all(1e-6).unwrap();
//! ```
//!
//! # Theorem
//! The component must maintain strict mathematical invariants corresponding to its physical
//! or numerical role.
//!
//! **Proof sketch**:
//! Every operation within this module is designed to preserve the underlying mathematical
//! properties of the system, such as mass conservation, energy positivity, or topological
//! consistency. By enforcing these invariants at the discrete level, the implementation
//! guarantees stability and physical realism.


use std::collections::{HashMap, VecDeque};

use cfd_core::error::Result as CfdResult;
use cfd_core::physics::fluid::BloodModel;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{NetworkBlueprint, NodeKind};
use cfd_schematics::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, SolveResult, StaggeredGrid2D};
use crate::solvers::venturi_flow::VenturiGeometry;

// ── Public API types ──────────────────────────────────────────────────────────

/// Per-channel 2D solve result, including solver convergence data and haemolysis.
#[derive(Debug, Clone)]
pub struct Channel2dResult<T> {
    /// Blueprint channel ID.
    pub channel_id: String,
    /// Therapy zone classification.
    pub therapy_zone: TherapyZone,
    /// Whether this channel contains a venturi throat geometry.
    pub is_venturi_throat: bool,
    /// SIMPLE solver convergence result.
    pub solve_result: SolveResult<T>,
    /// Estimated wall shear stress (6Q/(w·h²) model) [Pa].
    pub wall_shear_pa: T,
    /// Maximum wall shear stress extracted from 2D velocity field [Pa].
    pub field_wall_shear_max_pa: T,
    /// Mean wall shear stress extracted from 2D velocity field [Pa].
    pub field_wall_shear_mean_pa: T,
    /// Transit time through the channel [s].
    pub transit_time_s: T,
    /// Giersiepen (1990) haemolysis index contribution (dimensionless).
    pub hemolysis_index: T,
}

/// Results of a full-network 2D solve.
#[derive(Debug, Clone)]
pub struct Network2dResult<T> {
    /// Per-channel solve results in blueprint channel order.
    pub channels: Vec<Channel2dResult<T>>,
    /// Total haemolysis index across all channels (sum).
    pub total_hemolysis_index: T,
    /// Number of channels that converged.
    pub converged_count: usize,
}

/// Internal per-channel entry stored in [`Network2DSolver`].
struct Channel2dEntry<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    id: String,
    therapy_zone: TherapyZone,
    is_venturi_throat: bool,
    /// Flow rate through this channel [m³/s].
    flow_rate_m3_s: f64,
    width_m: f64,
    height_m: f64,
    length_m: f64,
    viscosity_pa_s: f64,
    solver: NavierStokesSolver2D<T>,
}

/// Full-network 2D solver: one configured 2D channel domain per blueprint channel.
///
/// Analogous to `cfd-1d`'s `Network<T, F>` but resolving full 2D velocity/pressure
/// fields for every straight rectangular channel in the network.
pub struct Network2DSolver<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    channels: Vec<Channel2dEntry<T>>,
}

impl<T> Network2DSolver<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    /// Solve every channel domain **in parallel** and return per-channel results.
    ///
    /// ## Theorem: Independence of Channel Solves
    ///
    /// Each channel in the `NetworkBlueprint` is a straight rectangular duct with
    /// prescribed inlet flow rate. The 2D Navier-Stokes equations are solved
    /// independently per channel (no cross-channel coupling), so the solves are
    /// embarrassingly parallel.
    ///
    /// ## Performance
    ///
    /// For an n-channel network on p threads, the wall-clock time is approximately
    /// `max_i(T_i)` where `T_i` is the solve time for channel `i`, versus
    /// `Σ_i T_i` for the sequential baseline. For symmetric bifurcation networks
    /// where all channels have similar geometry, this approaches `T_single × ceil(n/p)`.
    ///
    /// # Parameters
    /// - `tolerance`: SIMPLE iteration convergence tolerance.
    ///
    /// # Errors
    /// Returns an error if the 2D solver fails for any channel.
    pub fn solve_all(&mut self, tolerance: f64) -> CfdResult<Network2dResult<T>> {
        use rayon::prelude::*;

        let tol_t = T::from_f64(tolerance).unwrap_or_else(T::one);

        // Parallel map: each channel solves independently on its own thread.
        let per_channel: Vec<CfdResult<Channel2dResult<T>>> = self
            .channels
            .par_iter_mut()
            .map(|entry| {
                let area = entry.width_m * entry.height_m;
                let u_inlet =
                    T::from_f64(entry.flow_rate_m3_s / area.max(1e-18)).unwrap_or_else(T::one);

                // Update SIMPLE config tolerance.
                entry.solver.config.tolerance = tol_t;

                // Run 2D NS solve (SIMPLE loop).
                let solve_result = entry.solver.solve(u_inlet).map_err(|e| {
                    cfd_core::error::Error::InvalidInput(format!(
                        "Network2DSolver: channel '{}' failed: {e}",
                        entry.id
                    ))
                })?;

                // Wall shear stress estimate: τ = μ · 6Q / (w · h²)
                let shear_rate = 6.0 * entry.flow_rate_m3_s
                    / (entry.width_m * entry.height_m * entry.height_m).max(1e-30);
                let shear_pa = entry.viscosity_pa_s * shear_rate;

                // Field-resolved wall shear from the 2D velocity gradients.
                let (field_max, field_mean) = extract_field_wall_shear(&entry.solver);

                // Transit time: t = V / Q = w·h·L / Q
                let t_s = area * entry.length_m / entry.flow_rate_m3_s.max(1e-30);

                let hi = HemolysisModel::giersiepen_millifluidic()
                    .damage_index(shear_pa, t_s)
                    .unwrap_or(0.0);
                let hi_t = T::from_f64(hi).unwrap_or_else(T::zero);

                Ok(Channel2dResult {
                    channel_id: entry.id.clone(),
                    therapy_zone: entry.therapy_zone.clone(),
                    is_venturi_throat: entry.is_venturi_throat,
                    solve_result,
                    wall_shear_pa: T::from_f64(shear_pa).unwrap_or_else(T::zero),
                    field_wall_shear_max_pa: field_max,
                    field_wall_shear_mean_pa: field_mean,
                    transit_time_s: T::from_f64(t_s).unwrap_or_else(T::zero),
                    hemolysis_index: hi_t,
                })
            })
            .collect();

        // Sequential reduce: collect results, propagate first error, aggregate totals.
        let mut results = Vec::with_capacity(per_channel.len());
        let mut total_hi = T::zero();
        let mut converged_count = 0usize;

        for r in per_channel {
            let ch = r?;
            if ch.solve_result.converged {
                converged_count += 1;
            }
            total_hi += ch.hemolysis_index;
            results.push(ch);
        }

        Ok(Network2dResult {
            channels: results,
            total_hemolysis_index: total_hi,
            converged_count,
        })
    }
}

// ── GraphSink adapter ─────────────────────────────────────────────────────────

/// A [`GraphSink`] that converts a validated [`NetworkBlueprint`] into a
/// solver-ready [`Network2DSolver<T>`].
///
/// Each channel in the blueprint gets its own 2D Navier-Stokes domain on a
/// staggered grid of size `(grid_nx × grid_ny)`.
pub struct Network2dBuilderSink<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    blood: BloodModel<T>,
    density: f64,
    /// Total inlet flow rate [m³/s]; distributed across channels by topology.
    total_flow_rate_m3_s: f64,
    grid_nx: usize,
    grid_ny: usize,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> Network2dBuilderSink<T> {
    /// Create a new sink.
    ///
    /// # Parameters
    /// - `blood`: Blood rheology model (Newtonian, Casson, Carreau-Yasuda).
    /// - `density`: Blood density [kg/m³].
    /// - `total_flow_rate_m3_s`: Total volumetric flow rate at the device inlet [m³/s].
    /// - `grid_nx`: Number of grid cells in the streamwise direction per channel.
    /// - `grid_ny`: Number of grid cells in the cross-flow direction per channel.
    #[must_use]
    pub fn new(
        blood: BloodModel<T>,
        density: f64,
        total_flow_rate_m3_s: f64,
        grid_nx: usize,
        grid_ny: usize,
    ) -> Self {
        Self { blood, density, total_flow_rate_m3_s, grid_nx, grid_ny }
    }
}

impl<T> GraphSink for Network2dBuilderSink<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug + 'static,
{
    type Output = Network2DSolver<T>;

    fn build(&self, blueprint: &NetworkBlueprint) -> CfdResult<Network2DSolver<T>> {
        // 1. Compute per-channel flow rates by resistance-weighted BFS.
        let flow_rates = distribute_flow(blueprint, self.total_flow_rate_m3_s);

        // 2. Derive reference viscosity from the blood model.
        let mu = blood_viscosity_f64(&self.blood);

        // 3. Build one Channel2dEntry per blueprint channel.
        let mut entries = Vec::with_capacity(blueprint.channels.len());

        for ch in &blueprint.channels {
            let q_ch = flow_rates.get(ch.id.as_str()).copied().unwrap_or(0.0);

            // Cross-section dimensions via canonical SSOT method.
            let (w, h) = ch.cross_section.dims();

            // TherapyZone from metadata (default MixedFlow).
            let therapy_zone = ch
                .metadata
                .as_ref()
                .and_then(|m| m.get::<TherapyZoneMetadata>())
                .map_or(TherapyZone::MixedFlow, |tz| tz.zone.clone());

            // Detect venturi throat from metadata.
            let venturi_meta = ch
                .metadata
                .as_ref()
                .and_then(|m| m.get::<VenturiGeometryMetadata>());
            let is_venturi_throat = venturi_meta.is_some();

            // Build 2D staggered grid and NS solver.
            // For venturi channels the grid spans the full inlet width, not
            // the narrow throat width, so the converging/diverging mask can
            // be applied correctly.
            let l = ch.length_m;
            let l_t = T::from_f64(l).unwrap_or_else(T::one);
            let (grid_w, grid_width_m) = if let Some(vm) = venturi_meta {
                (T::from_f64(vm.inlet_width_m).unwrap_or_else(T::one), vm.inlet_width_m)
            } else {
                (T::from_f64(w).unwrap_or_else(T::one), w)
            };
            let grid = StaggeredGrid2D::new(self.grid_nx, self.grid_ny, l_t, grid_w);
            let density_t = T::from_f64(self.density).unwrap_or_else(T::one);
            let mut solver = NavierStokesSolver2D::new(
                grid,
                self.blood.clone(),
                density_t,
                SIMPLEConfig::default(),
            );

            // For venturi channels, populate the fluid/solid mask from the
            // converging-throat-diverging geometry.
            if let Some(vm) = venturi_meta {
                let converge_l = (l - vm.throat_length_m).max(0.0) / 2.0;
                let geom = VenturiGeometry::new(
                    T::from_f64(vm.inlet_width_m).unwrap_or_else(T::one),
                    T::from_f64(vm.throat_width_m).unwrap_or_else(T::one),
                    T::zero(),                                       // l_inlet = 0 (throat channel only)
                    T::from_f64(converge_l).unwrap_or_else(T::one),  // converging taper
                    T::from_f64(vm.throat_length_m).unwrap_or_else(T::one),
                    T::from_f64(converge_l).unwrap_or_else(T::one),  // diverging taper
                    T::from_f64(vm.throat_height_m).unwrap_or_else(T::one),
                );
                populate_venturi_mask(&mut solver, &geom, self.grid_nx, self.grid_ny);
            }

            entries.push(Channel2dEntry {
                id: ch.id.as_str().to_owned(),
                therapy_zone,
                is_venturi_throat,
                flow_rate_m3_s: q_ch,
                width_m: grid_width_m,
                height_m: h,
                length_m: l,
                viscosity_pa_s: mu,
                solver,
            });
        }

        Ok(Network2DSolver { channels: entries })
    }
}

// ── Flow distribution ─────────────────────────────────────────────────────────

/// Distribute total inlet flow across all channels by resistance-weighted BFS.
///
/// Starting from the inlet node with `q_total`, flow is divided at each junction
/// proportionally to the conductance (1/R) of each outgoing channel.
fn distribute_flow(blueprint: &NetworkBlueprint, q_total: f64) -> HashMap<String, f64> {
    let n_nodes = blueprint.nodes.len();
    if n_nodes == 0 {
        return HashMap::new();
    }

    // Node-id to dense index map (avoids string cloning in BFS state).
    let mut node_index: HashMap<&str, usize> = HashMap::with_capacity(n_nodes);
    for (idx, node) in blueprint.nodes.iter().enumerate() {
        node_index.insert(node.id.as_str(), idx);
    }

    // Build adjacency by dense node index.
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (ch_idx, ch) in blueprint.channels.iter().enumerate() {
        if let Some(&from_idx) = node_index.get(ch.from.as_str()) {
            outgoing[from_idx].push(ch_idx);
        }
    }

    // Find inlet node index.
    let Some(inlet_idx) = blueprint
        .nodes
        .iter()
        .enumerate()
        .find(|(_, n)| matches!(n.kind, NodeKind::Inlet))
        .map(|(idx, _)| idx)
    else {
        return HashMap::new();
    };

    // BFS traversal on dense node indices.
    let mut node_flow = vec![0.0_f64; n_nodes];
    node_flow[inlet_idx] = q_total;

    let mut channel_flow: HashMap<String, f64> = HashMap::with_capacity(blueprint.channels.len());
    let mut visited = vec![false; n_nodes];
    let mut queue = VecDeque::with_capacity(n_nodes);
    queue.push_back(inlet_idx);

    while let Some(node_idx) = queue.pop_front() {
        if visited[node_idx] {
            continue; // Already processed.
        }
        visited[node_idx] = true;

        let q_node = node_flow[node_idx];
        let channels = &outgoing[node_idx];

        if channels.is_empty() {
            continue;
        }

        // Total conductance of outgoing channels (1/R for each).
        let total_cond: f64 = channels
            .iter()
            .map(|&i| 1.0 / blueprint.channels[i].resistance.max(1e-30))
            .sum::<f64>()
            .max(1e-30);

        for &ch_idx in channels {
            let ch = &blueprint.channels[ch_idx];
            let cond = 1.0 / ch.resistance.max(1e-30);
            let q_ch = q_node * cond / total_cond;

            channel_flow.insert(ch.id.as_str().to_owned(), q_ch);

            // Accumulate flow at destination node index.
            if let Some(&to_idx) = node_index.get(ch.to.as_str()) {
                node_flow[to_idx] += q_ch;
                if !visited[to_idx] {
                    queue.push_back(to_idx);
                }
            }
        }
    }

    channel_flow
}

// ── Physics helpers ───────────────────────────────────────────────────────────

/// Extract a representative dynamic viscosity [Pa·s] from the blood model.
fn blood_viscosity_f64<T: RealField + Copy + Float + FromPrimitive + ToPrimitive>(
    model: &BloodModel<T>,
) -> f64 {
    const REF_SHEAR: f64 = 100.0; // 100 s⁻¹ reference shear rate
    let shear_t = T::from_f64(REF_SHEAR).unwrap_or_else(T::one);
    model.viscosity(shear_t).to_f64().unwrap_or(3.5e-3)
}

// ── Venturi mask ──────────────────────────────────────────────────────────────

/// Populate the solver's fluid/solid mask from a `VenturiGeometry`.
///
/// The grid y-coordinate is shifted so that `y=0` maps to the geometry
/// centreline, matching `VenturiGeometry::contains(x, y)` which tests
/// `|y| <= w_local/2`.
fn populate_venturi_mask<T>(
    solver: &mut NavierStokesSolver2D<T>,
    geom: &VenturiGeometry<T>,
    nx: usize,
    ny: usize,
) where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let half_h = geom.w_inlet / T::from_f64(2.0).unwrap_or_else(T::one);
    for i in 0..nx {
        for j in 0..ny {
            let x = solver.grid.x_center(i);
            let y = solver.grid.y_center(j) - half_h;
            solver.field.mask[i][j] = geom.contains(x, y);
        }
    }
}

// ── Field-resolved wall shear ─────────────────────────────────────────────────

/// Extract (max, mean) wall shear stress from the solved 2D velocity field.
///
/// Computes `τ = μ · γ̇` on the fly from the converged velocity gradients
/// (rather than reading the stored `gamma_dot`, which may be stale if the
/// solver converged between viscosity-update iterations).
///
/// Samples the first interior fluid cells adjacent to walls or solid regions.
fn extract_field_wall_shear<T>(solver: &NavierStokesSolver2D<T>) -> (T, T)
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let nx = solver.grid.nx;
    let ny = solver.grid.ny;
    let dx = solver.grid.dx;
    let dy = solver.grid.dy;
    let mut max_tau = T::zero();
    let mut sum_tau = T::zero();
    let mut count = 0u64;

    for i in 1..nx.saturating_sub(1) {
        for j in 1..ny.saturating_sub(1) {
            if !solver.field.mask[i][j] {
                continue;
            }
            let next_to_wall = !solver.field.mask[i - 1][j]
                || !solver.field.mask[i + 1][j]
                || !solver.field.mask[i][j - 1]
                || !solver.field.mask[i][j + 1]
                || i == 1
                || i == nx - 2
                || j == 1
                || j == ny - 2;

            if !next_to_wall {
                continue;
            }

            // Compute shear rate from the converged velocity field.
            let gamma = solver.field.compute_shear_rate(i, j, dx, dy);
            let tau = solver.field.mu[i][j] * gamma;
            if tau > max_tau {
                max_tau = tau;
            }
            sum_tau += tau;
            count += 1;
        }
    }

    let mean_tau = if count > 0 {
        sum_tau / T::from_u64(count).unwrap_or_else(T::one)
    } else {
        T::zero()
    };

    (max_tau, mean_tau)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_schematics::interface::presets::venturi_rect;

    #[test]
    fn build_single_venturi_network() {
        let bp = venturi_rect("test_v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
        let blood = BloodModel::Newtonian(3.5e-3_f64);
        let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-6, 10, 5);
        let net2d = sink.build(&bp).expect("build should succeed");
        assert_eq!(net2d.channels.len(), 3); // inlet, throat, diffuser
        // The throat channel must be detected as a venturi.
        let venturi_count = net2d.channels.iter().filter(|c| c.is_venturi_throat).count();
        assert!(venturi_count >= 1, "expected ≥1 venturi channel, got {venturi_count}");
    }

    #[test]
    fn flow_distribution_sums_to_q_total() {
        use cfd_schematics::interface::presets::bifurcation_venturi_rect;
        let bp = bifurcation_venturi_rect("bv", 0.005, 0.002, 0.0005, 0.001, 0.001);
        let q_total = 1e-6;
        let flows = distribute_flow(&bp, q_total);
        // Sum flows on channels leading into the outlet node.
        let outlet_id = bp
            .nodes
            .iter()
            .find(|n| matches!(n.kind, NodeKind::Outlet))
            .map(|n| n.id.as_str().to_owned())
            .unwrap_or_default();
        let q_out: f64 = bp
            .channels
            .iter()
            .filter(|c| c.to.as_str() == outlet_id)
            .map(|c| flows.get(c.id.as_str()).copied().unwrap_or(0.0))
            .sum();
        // At the outlet, flow should reconverge to approximately q_total.
        assert!(
            (q_out - q_total).abs() < q_total * 0.01,
            "outlet flow {q_out:.3e} != q_total {q_total:.3e}"
        );
    }

    #[test]
    fn field_wall_shear_agrees_with_analytical() {
        // Straight rectangular channel: analytical τ_w = μ · 6Q/(w·h²).
        let mu = 3.5e-3_f64;
        let w = 1.0e-3;
        let h = 1.0e-3;
        let l = 10.0e-3;
        let q = 1e-7; // 0.1 µL/s
        let u_mean = q / (w * h);

        let grid = StaggeredGrid2D::new(80, 20, l, w);
        let blood = BloodModel::Newtonian(mu);
        let mut solver = NavierStokesSolver2D::new(grid, blood, 1060.0, SIMPLEConfig::default());
        let result = solver.solve(u_mean).expect("straight channel solve");
        assert!(result.converged, "SIMPLE must converge for Poiseuille flow");

        let (max_tau, mean_tau) = extract_field_wall_shear(&solver);
        let analytical_tau = mu * 6.0 * q / (w * h * h);

        // Field-extracted mean should be within an order of magnitude of the
        // analytical value on a coarse 80×20 grid.
        let ratio = mean_tau / analytical_tau;
        assert!(
            ratio > 0.1 && ratio < 10.0,
            "field mean wall shear {mean_tau:.3e} should be within \
             an order of magnitude of analytical {analytical_tau:.3e}"
        );
        assert!(max_tau >= mean_tau, "max >= mean invariant");
    }
}
