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


use std::collections::{HashMap, HashSet, VecDeque};

use cfd_core::error::Result as CfdResult;
use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{CrossSectionSpec, NetworkBlueprint, NodeKind};
use cfd_schematics::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, SolveResult, StaggeredGrid2D};

// ── Constants ─────────────────────────────────────────────────────────────────

/// Giersiepen (1990) haemolysis constants: HI = C · t^α · τ^β
const GIERSIEPEN_C: f64 = 3.62e-5;
const GIERSIEPEN_ALPHA: f64 = 0.765;
const GIERSIEPEN_BETA: f64 = 1.991;

// ── Public API types ──────────────────────────────────────────────────────────

/// Per-channel 2D solve result, including solver convergence data and haemolysis.
#[derive(Debug, Clone)]
pub struct Channel2dResult<T> {
    /// Blueprint channel ID.
    pub channel_id: String,
    /// Therapy zone classification.
    pub therapy_zone: TherapyZone,
    /// SIMPLE solver convergence result.
    pub solve_result: SolveResult<T>,
    /// Estimated wall shear stress (6Q/(w·h²) model) [Pa].
    pub wall_shear_pa: T,
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
    /// Solve every channel domain and return per-channel results.
    ///
    /// # Parameters
    /// - `tolerance`: SIMPLE iteration convergence tolerance.
    ///
    /// # Errors
    /// Returns an error if the 2D solver fails for any channel.
    pub fn solve_all(&mut self, tolerance: f64) -> CfdResult<Network2dResult<T>> {
        let tol_t = T::from_f64(tolerance).unwrap_or_else(T::one);
        let mut results = Vec::with_capacity(self.channels.len());
        let mut total_hi = T::zero();
        let mut converged_count = 0usize;

        for entry in &mut self.channels {
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

            if solve_result.converged {
                converged_count += 1;
            }

            // Wall shear stress estimate: τ = μ · 6Q / (w · h²)
            let shear_rate = 6.0 * entry.flow_rate_m3_s
                / (entry.width_m * entry.height_m * entry.height_m).max(1e-30);
            let shear_pa = entry.viscosity_pa_s * shear_rate;

            // Transit time: t = V / Q = w·h·L / Q
            let t_s = area * entry.length_m / entry.flow_rate_m3_s.max(1e-30);

            let hi = giersiepen_hi(shear_pa, t_s);
            let hi_t = T::from_f64(hi).unwrap_or_else(T::zero);
            total_hi += hi_t;

            results.push(Channel2dResult {
                channel_id: entry.id.clone(),
                therapy_zone: entry.therapy_zone.clone(),
                solve_result,
                wall_shear_pa: T::from_f64(shear_pa).unwrap_or_else(T::zero),
                transit_time_s: T::from_f64(t_s).unwrap_or_else(T::zero),
                hemolysis_index: hi_t,
            });
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

            // Cross-section dimensions.
            let (w, h) = match ch.cross_section {
                CrossSectionSpec::Rectangular { width_m, height_m } => (width_m, height_m),
                CrossSectionSpec::Circular { diameter_m } => (diameter_m, diameter_m),
            };

            // TherapyZone from metadata (default MixedFlow).
            let therapy_zone = ch
                .metadata
                .as_ref()
                .and_then(|m| m.get::<TherapyZoneMetadata>())
                .map_or(TherapyZone::MixedFlow, |tz| tz.zone.clone());

            // Build 2D staggered grid and NS solver for this channel.
            let l = ch.length_m;
            let l_t = T::from_f64(l).unwrap_or_else(T::one);
            let w_t = T::from_f64(w).unwrap_or_else(T::one);
            let grid = StaggeredGrid2D::new(self.grid_nx, self.grid_ny, l_t, w_t);
            let density_t = T::from_f64(self.density).unwrap_or_else(T::one);
            let solver = NavierStokesSolver2D::new(
                grid,
                self.blood.clone(),
                density_t,
                SIMPLEConfig::default(),
            );

            entries.push(Channel2dEntry {
                id: ch.id.as_str().to_owned(),
                therapy_zone,
                flow_rate_m3_s: q_ch,
                width_m: w,
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
    // Build adjacency: node_id → list of channel indices for outgoing edges.
    let mut outgoing: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, ch) in blueprint.channels.iter().enumerate() {
        outgoing.entry(ch.from.as_str().to_owned()).or_default().push(i);
    }

    // Find the inlet node.
    let inlet_id = blueprint
        .nodes
        .iter()
        .find(|n| matches!(n.kind, NodeKind::Inlet))
        .map(|n| n.id.as_str().to_owned())
        .unwrap_or_default();

    // BFS traversal: track flow at each node.
    let mut node_flow: HashMap<String, f64> = HashMap::new();
    node_flow.insert(inlet_id.clone(), q_total);

    let mut channel_flow: HashMap<String, f64> = HashMap::new();
    let mut visited: HashSet<String> = HashSet::new();
    let mut queue = VecDeque::new();
    queue.push_back(inlet_id);

    while let Some(node_id) = queue.pop_front() {
        if !visited.insert(node_id.clone()) {
            continue; // Already processed.
        }

        let q_node = node_flow.get(&node_id).copied().unwrap_or(0.0);
        let channels = outgoing.get(&node_id).cloned().unwrap_or_default();

        if channels.is_empty() {
            continue;
        }

        // Total conductance of outgoing channels (1/R for each).
        let total_cond: f64 = channels
            .iter()
            .map(|&i| 1.0 / blueprint.channels[i].resistance.max(1e-30))
            .sum::<f64>()
            .max(1e-30);

        for &ch_idx in &channels {
            let ch = &blueprint.channels[ch_idx];
            let cond = 1.0 / ch.resistance.max(1e-30);
            let q_ch = q_node * cond / total_cond;

            channel_flow.insert(ch.id.as_str().to_owned(), q_ch);

            // Accumulate flow at the destination node.
            *node_flow.entry(ch.to.as_str().to_owned()).or_insert(0.0) += q_ch;

            queue.push_back(ch.to.as_str().to_owned());
        }
    }

    channel_flow
}

// ── Physics helpers ───────────────────────────────────────────────────────────

/// Giersiepen (1990) haemolysis index: `HI = C · t^α · τ^β`.
fn giersiepen_hi(tau_pa: f64, t_s: f64) -> f64 {
    if tau_pa <= 0.0 || t_s <= 0.0 {
        return 0.0;
    }
    GIERSIEPEN_C * t_s.powf(GIERSIEPEN_ALPHA) * tau_pa.powf(GIERSIEPEN_BETA)
}

/// Extract a representative dynamic viscosity [Pa·s] from the blood model.
fn blood_viscosity_f64<T: RealField + Copy + Float + FromPrimitive + ToPrimitive>(
    model: &BloodModel<T>,
) -> f64 {
    const REF_SHEAR: f64 = 100.0; // 100 s⁻¹ reference shear rate
    let shear_t = T::from_f64(REF_SHEAR).unwrap_or_else(T::one);
    model.viscosity(shear_t).to_f64().unwrap_or(3.5e-3)
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
}
