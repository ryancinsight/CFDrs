use std::collections::HashMap;

use cfd_1d::domain::network::{apply_blueprint_boundary_conditions, network_from_blueprint};
use cfd_1d::{NetworkSolver, SolverConfig};
use cfd_core::error::{Error, Result as CfdResult};
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_math::nonlinear_solver::{AndersonAccelerator, AndersonConfig, AndersonMethod};
use cfd_schematics::domain::model::{NetworkBlueprint, NodeKind};
use cfd_schematics::geometry::metadata::{
    BranchBoundaryMetadata, BranchBoundarySpecification, JunctionFamily, JunctionGeometryMetadata,
};
use nalgebra::{DVector, RealField};
use num_traits::{Float, FromPrimitive, ToPrimitive};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use rayon::prelude::*;

use super::channel::solve_channel_entry;
use super::reference::{
    build_reference_trace_from_solved_network, reference_fluid, NetworkReferenceTrace,
};
use super::types::{Channel2dEntry, Channel2dResult, CoupledNetwork2dResult, Network2dResult};
use super::Network2DSolver;
use cfd_1d::PrimarySolveDiagnostics;

const COUPLING_ANDERSON_DEPTH: usize = 4;
const MAX_COUPLING_ITERATIONS: usize = 24;
const MIN_LINEAR_RESISTANCE: f64 = 1e-12;
const COUPLING_CONVERGENCE_FLOOR_RELATIVE: f64 = 6e-4;
const COUPLING_WEIGHT_MIN: f64 = 0.75;
const COUPLING_WEIGHT_MAX: f64 = 1.35;
const COUPLING_AITKEN_RELAXATION_MIN: f64 = 0.05;
const COUPLING_AITKEN_RELAXATION_MAX: f64 = 1.5;
const COUPLING_AITKEN_DROP_TOLERANCE: f64 = 1e-12;

struct CoupledPassOutcome<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    network: cfd_1d::domain::network::Network<T, ConstantPropertyFluid<T>>,
    reference_trace: NetworkReferenceTrace<T>,
    channel_results: Vec<Channel2dResult<T>>,
    all_channels_converged: bool,
    max_relative_flow_change_pct: T,
    max_relative_resistance_change_pct: T,
}

struct ResistanceUpdateCandidate<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    edge_idx: EdgeIndex,
    current_linear: T,
}

struct CoupledResistanceMixer<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    anderson: AndersonAccelerator<T>,
    previous_residual: Option<DVector<T>>,
    previous_relaxation: Option<DVector<T>>,
    relaxation_floor: T,
    relaxation_ceiling: T,
    residual_drop_tolerance: T,
}

impl<T> CoupledResistanceMixer<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    fn new(history_depth: usize) -> Self {
        Self {
            anderson: AndersonAccelerator::new(AndersonConfig::<T> {
                history_depth,
                relaxation: T::one(),
                drop_tolerance: T::from_f64(COUPLING_AITKEN_DROP_TOLERANCE)
                    .unwrap_or_else(T::default_epsilon),
                method: AndersonMethod::QR,
            }),
            previous_residual: None,
            previous_relaxation: None,
            relaxation_floor: T::from_f64(COUPLING_AITKEN_RELAXATION_MIN)
                .expect("analytical constant conversion"),
            relaxation_ceiling: T::from_f64(COUPLING_AITKEN_RELAXATION_MAX)
                .expect("analytical constant conversion"),
            residual_drop_tolerance: T::from_f64(COUPLING_AITKEN_DROP_TOLERANCE)
                .unwrap_or_else(T::default_epsilon),
        }
    }

    fn mix(&mut self, current: &DVector<T>, target: &DVector<T>) -> DVector<T> {
        if current.len() != target.len() || !vector_is_finite(current) || !vector_is_finite(target)
        {
            self.reset();
            return target.clone();
        }

        let aitken_candidate = self.vector_aitken_candidate(current, target);
        if !vector_is_finite(&aitken_candidate) {
            self.reset();
            return target.clone();
        }

        let anderson_candidate = self.anderson.compute_next(current, &aitken_candidate);
        if vector_is_finite(&anderson_candidate) {
            anderson_candidate
        } else {
            self.anderson.reset();
            aitken_candidate
        }
    }

    fn reset(&mut self) {
        self.previous_residual = None;
        self.previous_relaxation = None;
        self.anderson.reset();
    }

    fn vector_aitken_candidate(&mut self, current: &DVector<T>, target: &DVector<T>) -> DVector<T> {
        let residual = target - current;
        let mut relaxation = DVector::from_element(residual.len(), T::one());

        if let (Some(previous_residual), Some(previous_relaxation)) =
            (&self.previous_residual, &self.previous_relaxation)
        {
            if previous_residual.len() == residual.len()
                && previous_relaxation.len() == residual.len()
            {
                let residual_delta = &residual - previous_residual;
                for i in 0..residual.len() {
                    let previous_factor = previous_relaxation[i];
                    let delta = residual_delta[i];
                    relaxation[i] = if delta.is_finite()
                        && Float::abs(delta) > self.residual_drop_tolerance
                        && previous_factor.is_finite()
                    {
                        let raw = -previous_factor * previous_residual[i] / delta;
                        clamp_relaxation(raw, self.relaxation_floor, self.relaxation_ceiling)
                    } else {
                        clamp_relaxation(
                            previous_factor,
                            self.relaxation_floor,
                            self.relaxation_ceiling,
                        )
                    };
                }
            } else {
                self.reset();
            }
        }

        self.previous_residual = Some(residual.clone());
        self.previous_relaxation = Some(relaxation.clone());
        current + residual.component_mul(&relaxation)
    }
}

impl<T> Network2DSolver<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    /// Solve the junction-aware network using coupled 1D/2D iterations.
    ///
    /// The outer loop uses a vector Anderson mixer on the full resistance state,
    /// with branch/run asymmetry derived from explicit junction metadata in the
    /// blueprint. That keeps the coupled pass stable while letting split and
    /// merge branches relax at different rates.
    pub fn solve_coupled(&mut self, tolerance: f64) -> CfdResult<CoupledNetwork2dResult<T>> {
        let target_total_flow_m3_s = self
            .reference_trace
            .total_inlet_flow_m3_s
            .to_f64()
            .unwrap_or(0.0);
        let fluid =
            reference_fluid::<T>(self.reference_density_kg_m3, self.reference_viscosity_pa_s)?;
        let mut working_network =
            network_from_blueprint::<T, _>(&self.blueprint, fluid).map_err(|e| {
                Error::InvalidInput(format!(
                    "Network2DSolver coupled solve failed to build the 1D network: {e}"
                ))
            })?;
        let node_indices: HashMap<String, NodeIndex> = working_network
            .graph
            .node_indices()
            .filter_map(|idx| {
                working_network
                    .graph
                    .node_weight(idx)
                    .map(|node| (node.id.clone(), idx))
            })
            .collect();
        apply_blueprint_boundary_conditions(
            &mut working_network,
            &self.blueprint,
            &node_indices,
            T::one(),
            T::zero(),
        )?;
        let seed_state_by_channel_id: HashMap<&str, (T, T)> = self
            .channels
            .iter()
            .map(|entry| {
                let seed_flow =
                    T::from_f64(entry.flow_rate_m3_s).expect("analytical constant conversion");
                let seed_reference_flow = entry.reference_trace.flow_rate_m3_s;
                let seed_reference_drop = entry.reference_trace.pressure_drop_pa;
                let seed_resistance = if Float::abs(seed_reference_flow)
                    > T::from_f64(MIN_LINEAR_RESISTANCE).expect("analytical constant conversion")
                {
                    Float::abs(seed_reference_drop / seed_reference_flow)
                } else {
                    entry.reference_trace.resistance_pa_s_per_m3
                };
                (entry.id.as_str(), (seed_flow, seed_resistance))
            })
            .collect();
        for edge_idx in working_network.graph.edge_indices().collect::<Vec<_>>() {
            let edge_id = working_network
                .graph
                .edge_weight(edge_idx)
                .map(|edge| edge.id.clone());
            if let Some(edge) = working_network.graph.edge_weight_mut(edge_idx) {
                if !edge.resistance.is_finite() || edge.resistance <= T::zero() {
                    edge.resistance =
                        T::from_f64(MIN_LINEAR_RESISTANCE).expect("analytical constant conversion");
                }
                if let Some(edge_id) = edge_id.as_deref() {
                    if let Some((seed_flow, seed_resistance)) =
                        seed_state_by_channel_id.get(edge_id)
                    {
                        working_network.flow_rates[edge_idx.index()] = *seed_flow;
                        edge.flow_rate = *seed_flow;
                        edge.resistance = *seed_resistance;
                    }
                }
                edge.quad_coeff = T::zero();
            }
            if let Some(props) = working_network.properties.get_mut(&edge_idx) {
                if !props.resistance.is_finite() || props.resistance <= T::zero() {
                    props.resistance =
                        T::from_f64(MIN_LINEAR_RESISTANCE).expect("analytical constant conversion");
                }
                if let Some(edge_id) = edge_id.as_deref() {
                    if let Some((_, seed_resistance)) = seed_state_by_channel_id.get(edge_id) {
                        props.resistance = *seed_resistance;
                    }
                }
                props.resistance_update_policy =
                    cfd_1d::domain::network::ResistanceUpdatePolicy::FlowInvariant;
            }
        }

        let edge_index_by_id = edge_index_by_id(&working_network);
        let solver = NetworkSolver::<T, ConstantPropertyFluid<T>>::with_config(SolverConfig::<T> {
            tolerance: T::from_f64(1e-12).unwrap_or_else(T::default_epsilon),
            max_iterations: 500,
        });
        let mut coupling_mixer = CoupledResistanceMixer::<T>::new(COUPLING_ANDERSON_DEPTH);
        let channel_coupling_weights = build_channel_coupling_weights::<T>(&self.blueprint);
        let convergence_threshold_pct =
            T::from_f64(tolerance.max(COUPLING_CONVERGENCE_FLOOR_RELATIVE) * 100.0)
                .expect("analytical constant conversion");
        let mut previous_reference_trace: Option<NetworkReferenceTrace<T>> = None;
        let mut final_reference_trace = self.reference_trace.clone();
        let mut final_channel_results = Vec::new();
        let mut coupling_iterations = 0usize;
        let mut converged = false;
        let mut max_relative_flow_change_pct = T::zero();
        let mut max_relative_resistance_change_pct = T::zero();

        for iteration in 0..MAX_COUPLING_ITERATIONS {
            let outcome = run_coupled_pass(
                &solver,
                &self.blueprint,
                &mut self.channels,
                &edge_index_by_id,
                working_network,
                tolerance,
                target_total_flow_m3_s,
                &channel_coupling_weights,
                &mut coupling_mixer,
                previous_reference_trace.as_ref(),
                &self.reference_trace,
                self.separation_tracking_enabled,
            )?;

            coupling_iterations = iteration + 1;
            max_relative_flow_change_pct = outcome.max_relative_flow_change_pct;
            max_relative_resistance_change_pct = outcome.max_relative_resistance_change_pct;
            final_reference_trace = outcome.reference_trace.clone();
            final_channel_results = outcome.channel_results;
            working_network = outcome.network;

            converged = outcome.all_channels_converged
                && max_relative_flow_change_pct <= convergence_threshold_pct
                && max_relative_resistance_change_pct <= convergence_threshold_pct;

            if converged {
                break;
            }

            previous_reference_trace = Some(final_reference_trace.clone());
        }

        self.reference_trace = final_reference_trace.clone();
        for (entry, channel_trace) in self
            .channels
            .iter_mut()
            .zip(final_reference_trace.channel_traces.iter())
        {
            entry.reference_trace = channel_trace.clone();
            entry.flow_rate_m3_s = channel_trace.flow_rate_m3_s.to_f64().unwrap_or(0.0);
        }

        let mut results = Vec::with_capacity(final_channel_results.len());
        let mut total_hi = T::zero();
        let mut total_outlet_error_pct = T::zero();
        let mut max_outlet_error_pct = T::zero();
        let mut converged_count = 0usize;

        for channel in final_channel_results {
            if channel.solve_result.converged {
                converged_count += 1;
            }
            if channel.field_outlet_flow_error_pct > max_outlet_error_pct {
                max_outlet_error_pct = channel.field_outlet_flow_error_pct;
            }
            total_outlet_error_pct += channel.field_outlet_flow_error_pct;
            total_hi += channel.hemolysis_index;
            results.push(channel);
        }

        let mean_outlet_error_pct = if results.is_empty() {
            T::zero()
        } else {
            total_outlet_error_pct
                / T::from_usize(results.len()).expect("analytical constant conversion")
        };

        let result = Network2dResult {
            channels: results,
            total_hemolysis_index: total_hi,
            converged_count,
            max_field_outlet_flow_error_pct: max_outlet_error_pct,
            mean_field_outlet_flow_error_pct: mean_outlet_error_pct,
            reference_trace: final_reference_trace,
        };

        Ok(CoupledNetwork2dResult {
            result,
            projection: self.projection_summary_ref().clone(),
            coupling_iterations,
            converged,
            max_relative_flow_change_pct,
            max_relative_resistance_change_pct,
        })
    }
}

fn run_coupled_pass<T>(
    solver: &NetworkSolver<T, ConstantPropertyFluid<T>>,
    blueprint: &NetworkBlueprint,
    channels: &mut [Channel2dEntry<T>],
    edge_index_by_id: &HashMap<String, EdgeIndex>,
    working_network: cfd_1d::domain::network::Network<T, ConstantPropertyFluid<T>>,
    tolerance: f64,
    target_total_flow_m3_s: f64,
    channel_coupling_weights: &[T],
    coupling_mixer: &mut CoupledResistanceMixer<T>,
    previous_reference_trace: Option<&NetworkReferenceTrace<T>>,
    seed_reference_trace: &NetworkReferenceTrace<T>,
    separation_tracking_enabled: bool,
) -> CfdResult<CoupledPassOutcome<T>>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    let solve_input = working_network.clone();
    let solve_result = solver.solve_owned_network_with_diagnostics(solve_input);
    let (mut solved_network, diagnostics) = if let Ok(result) = solve_result {
        result
    } else {
        let fallback_network = working_network;
        (fallback_network, PrimarySolveDiagnostics::default())
    };
    let mut fell_back_to_seed_reference = false;
    let current_reference_trace = if let Ok(trace) = build_reference_trace_from_solved_network(
        blueprint,
        &solved_network,
        diagnostics,
        T::one(),
        target_total_flow_m3_s,
    ) {
        trace
    } else {
        fell_back_to_seed_reference = true;
        seed_reference_trace.clone()
    };

    let per_channel: Vec<CfdResult<Channel2dResult<T>>> = channels
        .par_iter_mut()
        .zip(current_reference_trace.channel_traces.par_iter())
        .map(|(entry, channel_trace)| {
            solve_channel_entry(
                entry,
                tolerance,
                channel_trace.flow_rate_m3_s.to_f64().unwrap_or(0.0),
                channel_trace,
                separation_tracking_enabled,
            )
        })
        .collect();

    let mut channel_results = Vec::with_capacity(per_channel.len());
    let mut all_channels_converged = true;
    for result in per_channel {
        let channel = result?;
        all_channels_converged &= channel.solve_result.converged;
        channel_results.push(channel);
    }

    if channel_coupling_weights.len() != channel_results.len() {
        return Err(Error::InvalidInput(format!(
            "Network2DSolver coupled solve expected {} channel coupling weights but received {}",
            channel_results.len(),
            channel_coupling_weights.len()
        )));
    }

    let flow_floor = T::from_f64(1e-30).expect("analytical constant conversion");
    let min_linear_resistance =
        T::from_f64(MIN_LINEAR_RESISTANCE).expect("analytical constant conversion");
    let mut current_linear_values = Vec::with_capacity(channel_results.len());
    let mut weighted_target_values = Vec::with_capacity(channel_results.len());
    let mut update_candidates = Vec::with_capacity(channel_results.len());
    let mut max_relative_resistance_change = T::zero();

    for ((channel_result, channel_trace), weight) in channel_results
        .iter()
        .zip(current_reference_trace.channel_traces.iter())
        .zip(channel_coupling_weights.iter())
    {
        let Some(&edge_idx) = edge_index_by_id.get(channel_result.channel_id.as_str()) else {
            return Err(Error::InvalidInput(format!(
                "Network2DSolver coupled solve missing edge '{}' in the 1D network",
                channel_result.channel_id.as_str()
            )));
        };

        let (current_linear, current_quad) = if fell_back_to_seed_reference {
            (channel_trace.resistance_pa_s_per_m3, T::zero())
        } else {
            let Some(edge_weight) = solved_network.graph.edge_weight(edge_idx) else {
                return Err(Error::InvalidInput(format!(
                    "Network2DSolver coupled solve missing edge weight for '{}'",
                    channel_result.channel_id.as_str()
                )));
            };
            (edge_weight.resistance, edge_weight.quad_coeff)
        };
        let flow_abs = Float::abs(channel_trace.flow_rate_m3_s);
        let target_linear = if flow_abs > flow_floor {
            let total_effective = channel_result.field_effective_resistance_pa_s_per_m3;
            Float::max(
                total_effective - current_quad * flow_abs,
                min_linear_resistance,
            )
        } else {
            current_linear
        };
        let weighted_target_linear = current_linear + (*weight * (target_linear - current_linear));

        current_linear_values.push(current_linear);
        weighted_target_values.push(weighted_target_linear);
        update_candidates.push(ResistanceUpdateCandidate {
            edge_idx,
            current_linear,
        });
    }

    let current_linear = DVector::from_vec(current_linear_values);
    let mut weighted_target = DVector::from_vec(weighted_target_values);
    clamp_vector_to_floor(&mut weighted_target, min_linear_resistance);
    if !vector_is_finite(&current_linear) || !vector_is_finite(&weighted_target) {
        return Err(Error::InvalidInput(
            "Network2DSolver coupled solve produced a non-finite resistance state".to_string(),
        ));
    }

    let mut next_linear = coupling_mixer.mix(&current_linear, &weighted_target);
    clamp_vector_to_floor(&mut next_linear, min_linear_resistance);
    if !vector_is_finite(&next_linear) {
        coupling_mixer.reset();
        next_linear = weighted_target.clone();
        clamp_vector_to_floor(&mut next_linear, min_linear_resistance);
    }

    for ((candidate, channel_trace), next_linear_value) in update_candidates
        .iter()
        .zip(current_reference_trace.channel_traces.iter())
        .zip(next_linear.iter())
    {
        let denom = Float::max(
            Float::max(
                Float::abs(candidate.current_linear),
                Float::abs(*next_linear_value),
            ),
            min_linear_resistance,
        );
        let relative_change = Float::abs(*next_linear_value - candidate.current_linear) / denom;
        if relative_change > max_relative_resistance_change {
            max_relative_resistance_change = relative_change;
        }

        if let Some(edge) = solved_network.graph.edge_weight_mut(candidate.edge_idx) {
            edge.resistance = *next_linear_value;
            edge.flow_rate = channel_trace.flow_rate_m3_s;
        }
        if let Some(props) = solved_network.properties.get_mut(&candidate.edge_idx) {
            props.resistance = *next_linear_value;
            props.resistance_update_policy =
                cfd_1d::domain::network::ResistanceUpdatePolicy::FlowInvariant;
        }
    }

    let mut max_relative_flow_change = T::zero();
    if let Some(previous) = previous_reference_trace {
        for (current, previous) in current_reference_trace
            .channel_traces
            .iter()
            .zip(previous.channel_traces.iter())
        {
            let denom = Float::max(
                Float::max(
                    Float::abs(current.flow_rate_m3_s),
                    Float::abs(previous.flow_rate_m3_s),
                ),
                flow_floor,
            );
            let change = Float::abs(current.flow_rate_m3_s - previous.flow_rate_m3_s) / denom;
            if change > max_relative_flow_change {
                max_relative_flow_change = change;
            }
        }
    }

    Ok(CoupledPassOutcome {
        network: solved_network,
        reference_trace: current_reference_trace,
        channel_results,
        all_channels_converged,
        max_relative_flow_change_pct: max_relative_flow_change
            * T::from_f64(100.0).expect("analytical constant conversion"),
        max_relative_resistance_change_pct: max_relative_resistance_change
            * T::from_f64(100.0).expect("analytical constant conversion"),
    })
}

pub(crate) fn build_channel_coupling_weights<T>(blueprint: &NetworkBlueprint) -> Vec<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    let node_lookup: HashMap<
        &str,
        (
            &NodeKind,
            Option<&JunctionGeometryMetadata>,
            Option<&BranchBoundaryMetadata>,
        ),
    > = blueprint
        .nodes
        .iter()
        .map(|node| {
            (
                node.id.as_str(),
                (
                    &node.kind,
                    node.metadata
                        .as_ref()
                        .and_then(|metadata| metadata.get::<JunctionGeometryMetadata>()),
                    node.metadata
                        .as_ref()
                        .and_then(|metadata| metadata.get::<BranchBoundaryMetadata>()),
                ),
            )
        })
        .collect();

    blueprint
        .channels
        .iter()
        .map(|channel| {
            let mut endpoint_weights = Vec::with_capacity(2);
            if let Some(weight) =
                endpoint_coupling_weight(&node_lookup, channel.from.as_str(), true)
            {
                endpoint_weights.push(weight);
            }
            if let Some(weight) = endpoint_coupling_weight(&node_lookup, channel.to.as_str(), false)
            {
                endpoint_weights.push(weight);
            }

            let combined_weight = if endpoint_weights.is_empty() {
                1.0
            } else {
                endpoint_weights.iter().copied().sum::<f64>() / endpoint_weights.len() as f64
            };
            T::from_f64(clamp_coupling_weight(combined_weight))
                .expect("analytical constant conversion")
        })
        .collect()
}

fn endpoint_coupling_weight(
    node_lookup: &HashMap<
        &str,
        (
            &NodeKind,
            Option<&JunctionGeometryMetadata>,
            Option<&BranchBoundaryMetadata>,
        ),
    >,
    node_id: &str,
    is_outgoing: bool,
) -> Option<f64> {
    let (node_kind, junction_meta, branch_boundary_meta) = *node_lookup.get(node_id)?;

    let geometry_weight = if *node_kind == NodeKind::Junction {
        junction_meta.map(|meta| {
            let branch_edge = match meta.junction_family {
                JunctionFamily::Merge => !is_outgoing,
                JunctionFamily::Bifurcation
                | JunctionFamily::Trifurcation
                | JunctionFamily::Tee => is_outgoing,
                JunctionFamily::Cross => true,
            };

            let angle_deg = if branch_edge {
                representative_abs_angle_deg(&meta.branch_angles_deg)
                    .or_else(|| representative_abs_angle_deg(&meta.merge_angles_deg))
            } else {
                representative_abs_angle_deg(&meta.merge_angles_deg)
                    .or_else(|| representative_abs_angle_deg(&meta.branch_angles_deg))
            }
            .unwrap_or(0.0);

            match meta.junction_family {
                JunctionFamily::Bifurcation => {
                    branch_run_weight(branch_edge, angle_deg, 1.18, 0.90, 0.12, 0.05)
                }
                JunctionFamily::Trifurcation => {
                    branch_run_weight(branch_edge, angle_deg, 1.15, 0.92, 0.10, 0.05)
                }
                JunctionFamily::Tee => {
                    branch_run_weight(branch_edge, angle_deg, 1.12, 0.94, 0.08, 0.04)
                }
                JunctionFamily::Merge => {
                    branch_run_weight(branch_edge, angle_deg, 1.20, 0.88, 0.12, 0.05)
                }
                JunctionFamily::Cross => 1.0,
            }
        })
    } else {
        None
    };

    let boundary_weight =
        branch_boundary_meta.map(|meta| boundary_coupling_weight(meta, is_outgoing));

    match (geometry_weight, boundary_weight) {
        (Some(geometry_weight), Some(boundary_weight)) => {
            Some(clamp_coupling_weight(geometry_weight * boundary_weight))
        }
        (Some(geometry_weight), None) => Some(geometry_weight),
        (None, Some(boundary_weight)) => Some(boundary_weight),
        (None, None) => None,
    }
}

fn branch_run_weight(
    branch_edge: bool,
    angle_deg: f64,
    branch_base: f64,
    run_base: f64,
    branch_angle_gain: f64,
    run_angle_gain: f64,
) -> f64 {
    let angle_ratio = (angle_deg.abs() / 90.0).clamp(0.0, 1.0);
    let raw = if branch_edge {
        branch_base * (1.0 + branch_angle_gain * angle_ratio)
    } else {
        run_base * (1.0 - run_angle_gain * angle_ratio)
    };
    clamp_coupling_weight(raw)
}

fn boundary_coupling_weight(boundary: &BranchBoundaryMetadata, is_outgoing: bool) -> f64 {
    let (source_like, source_strength, sink_strength) = match boundary.boundary {
        BranchBoundarySpecification::Pressure { pressure_pa } if pressure_pa.is_finite() => {
            if pressure_pa == 0.0 {
                return 1.0;
            }
            (pressure_pa > 0.0, 1.04, 0.96)
        }
        BranchBoundarySpecification::FlowRate { flow_rate_m3_s } if flow_rate_m3_s.is_finite() => {
            if flow_rate_m3_s == 0.0 {
                return 1.0;
            }
            (flow_rate_m3_s > 0.0, 1.08, 0.92)
        }
        _ => return 1.0,
    };

    let raw = if source_like {
        if is_outgoing {
            source_strength
        } else {
            sink_strength
        }
    } else if is_outgoing {
        sink_strength
    } else {
        source_strength
    };
    clamp_coupling_weight(raw)
}

fn clamp_coupling_weight(value: f64) -> f64 {
    value.clamp(COUPLING_WEIGHT_MIN, COUPLING_WEIGHT_MAX)
}

fn representative_abs_angle_deg(angles: &[f64]) -> Option<f64> {
    angles
        .iter()
        .copied()
        .filter(|angle| angle.is_finite())
        .map(f64::abs)
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
}

fn vector_is_finite<T>(vector: &DVector<T>) -> bool
where
    T: Float,
{
    vector.iter().all(|value| value.is_finite())
}

fn clamp_relaxation<T>(value: T, floor: T, ceiling: T) -> T
where
    T: Float + Copy,
{
    if !value.is_finite() {
        return floor;
    }
    Float::max(Float::min(value, ceiling), floor)
}

fn clamp_vector_to_floor<T>(vector: &mut DVector<T>, floor: T)
where
    T: Float + Copy,
{
    for value in vector.iter_mut() {
        if !value.is_finite() || *value <= floor {
            *value = floor;
        }
    }
}

fn edge_index_by_id<T>(
    network: &cfd_1d::domain::network::Network<T, ConstantPropertyFluid<T>>,
) -> HashMap<String, EdgeIndex>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    network
        .graph
        .edge_references()
        .map(|edge_ref| (edge_ref.weight().id.clone(), edge_ref.id()))
        .collect()
}
