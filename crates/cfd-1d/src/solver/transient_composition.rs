//! Transient inlet-composition simulation for 1D networks.
//!
//! This module enables time-varying inlet fluid composition schedules and computes
//! node/edge mixture states over time using existing solved network flow rates.

use crate::network::Network;
use super::{NetworkProblem, NetworkSolver};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Mixture composition keyed by `fluid_id` with mass/volume fractions.
#[derive(Debug, Clone)]
pub struct MixtureComposition<T: RealField + Copy> {
    /// Fraction per fluid id.
    pub fractions: HashMap<i32, T>,
}

impl<T: RealField + Copy + FromPrimitive> MixtureComposition<T> {
    /// Create a new mixture and normalize to unit sum (when non-empty).
    #[must_use]
    pub fn new(mut fractions: HashMap<i32, T>) -> Self {
        let sum = fractions.values().copied().fold(T::zero(), |acc, v| acc + v);
        if sum > T::zero() {
            for value in fractions.values_mut() {
                *value /= sum;
            }
        }
        Self { fractions }
    }

    /// Empty mixture.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            fractions: HashMap::new(),
        }
    }

    /// Weighted blend of incoming mixtures.
    #[must_use]
    pub fn blend_weighted(inputs: &[(Self, T)]) -> Self {
        if inputs.is_empty() {
            return Self::empty();
        }

        let total_weight = inputs
            .iter()
            .map(|(_, w)| *w)
            .fold(T::zero(), |acc, v| acc + v);

        if total_weight <= T::zero() {
            return Self::empty();
        }

        let mut blended: HashMap<i32, T> = HashMap::new();
        for (mixture, weight) in inputs {
            for (fluid_id, frac) in &mixture.fractions {
                let contribution = (*frac) * (*weight) / total_weight;
                let entry = blended.entry(*fluid_id).or_insert(T::zero());
                *entry += contribution;
            }
        }

        Self::new(blended)
    }

    /// Compare with tolerance.
    #[must_use]
    pub fn approximately_equals(&self, other: &Self, tolerance: T) -> bool {
        let mut keys: Vec<i32> = self.fractions.keys().copied().collect();
        for key in other.fractions.keys() {
            if !keys.contains(key) {
                keys.push(*key);
            }
        }

        keys.into_iter().all(|k| {
            let a = *self.fractions.get(&k).unwrap_or(&T::zero());
            let b = *other.fractions.get(&k).unwrap_or(&T::zero());
            (a - b).abs() <= tolerance
        })
    }
}

/// Piecewise-constant inlet mixture event.
#[derive(Debug, Clone)]
pub struct InletCompositionEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Node index where the inlet composition applies.
    pub node_index: usize,
    /// Mixture after this event.
    pub mixture: MixtureComposition<T>,
}

/// Piecewise-constant edge flow event for transient pump-style control.
#[derive(Debug, Clone)]
pub struct EdgeFlowEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Edge index whose flow is updated.
    pub edge_index: usize,
    /// Flow rate after this event.
    pub flow_rate: T,
}

/// Piecewise-constant node pressure boundary event for transient pressure-pump control.
#[derive(Debug, Clone)]
pub struct PressureBoundaryEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Node index whose pressure boundary is updated.
    pub node_index: usize,
    /// Pressure value after this event.
    pub pressure: T,
}

/// Composition state at one simulation timepoint.
#[derive(Debug, Clone)]
pub struct CompositionState<T: RealField + Copy> {
    /// Simulation time.
    pub time: T,
    /// Node mixture compositions keyed by node index.
    pub node_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge mixture compositions keyed by edge index.
    pub edge_mixtures: HashMap<usize, MixtureComposition<T>>,
    /// Edge flow rates keyed by edge index for this snapshot.
    pub edge_flow_rates: HashMap<usize, T>,
}

impl<T: RealField + Copy> CompositionState<T> {
    /// Return average fluid concentrations in an edge at this state.
    ///
    /// In the current architecture, each edge stores a single mixed composition
    /// per snapshot, so this returns that composition as fluid concentrations.
    #[must_use]
    pub fn average_fluid_concentrations_in_edge(&self, edge_index: usize) -> Option<HashMap<i32, T>> {
        self.edge_mixtures
            .get(&edge_index)
            .map(|mixture| mixture.fractions.clone())
    }
}

/// Time control configuration for transient simulations.
///
/// This mirrors MMFT-style controls:
/// - `duration`: total simulation time
/// - `result_time_step`: interval between returned states
/// - `calculation_time_step`: internal update interval
#[derive(Debug, Clone)]
pub struct SimulationTimeConfig<T: RealField + Copy> {
    /// Total simulation duration.
    pub duration: T,
    /// Interval between result snapshots.
    pub result_time_step: T,
    /// Internal calculation/update timestep.
    pub calculation_time_step: T,
}

impl<T: RealField + Copy + FromPrimitive> SimulationTimeConfig<T> {
    /// Create a new time control configuration.
    #[must_use]
    pub fn new(duration: T, result_time_step: T, calculation_time_step: T) -> Self {
        Self {
            duration,
            result_time_step,
            calculation_time_step,
        }
    }

    fn validate(&self) -> Result<()> {
        if self.duration < T::zero() {
            return Err(Error::InvalidInput(
                "Simulation duration must be non-negative".to_string(),
            ));
        }
        if self.result_time_step <= T::zero() {
            return Err(Error::InvalidInput(
                "Simulation result timestep must be positive".to_string(),
            ));
        }
        if self.calculation_time_step <= T::zero() {
            return Err(Error::InvalidInput(
                "Simulation calculation timestep must be positive".to_string(),
            ));
        }
        Ok(())
    }

    /// Generate output snapshot timepoints in `[0, duration]`.
    pub fn result_timepoints(&self) -> Result<Vec<T>> {
        self.validate()?;
        Ok(Self::uniform_timepoints(self.duration, self.result_time_step))
    }

    /// Generate internal calculation timepoints in `[0, duration]`.
    pub fn calculation_timepoints(&self) -> Result<Vec<T>> {
        self.validate()?;
        Ok(Self::uniform_timepoints(
            self.duration,
            self.calculation_time_step,
        ))
    }

    fn uniform_timepoints(duration: T, dt: T) -> Vec<T> {
        let tolerance = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);
        if duration <= tolerance {
            return vec![T::zero()];
        }

        let mut points = vec![T::zero()];
        let mut t = dt;
        while t < duration {
            points.push(t);
            t += dt;
        }

        let needs_terminal = points
            .last()
            .map(|last| (*last - duration).abs() > tolerance)
            .unwrap_or(true);

        if needs_terminal {
            points.push(duration);
        }

        points
    }
}

/// Simulator for transient inlet composition switching with instantaneous node mixing.
pub struct TransientCompositionSimulator;

impl TransientCompositionSimulator {
    /// Simulate using MMFT-style timing controls.
    pub fn simulate_with_time_config<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletCompositionEvent<T>>,
        timing: SimulationTimeConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        let result_timepoints = timing.result_timepoints()?;
        let mut all_timepoints = result_timepoints.clone();
        all_timepoints.extend(timing.calculation_timepoints()?);

        for event in &events {
            if event.time >= T::zero() && event.time <= timing.duration {
                all_timepoints.push(event.time);
            }
        }

        let tolerance = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);
        let merged_timepoints = Self::sort_unique_timepoints(all_timepoints, tolerance);
        let all_states = Self::simulate(network, events, merged_timepoints)?;

        let mut sampled_states = Vec::with_capacity(result_timepoints.len());
        for target_time in result_timepoints {
            let state = all_states
                .iter()
                .find(|state| Self::times_close(state.time, target_time, tolerance))
                .cloned()
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "Failed to sample one or more configured result timepoints".to_string(),
                    )
                })?;
            sampled_states.push(state);
        }

        Ok(sampled_states)
    }

    /// Simulate composition over the provided timepoints using the network flow field.
    pub fn simulate<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_with_flow_events(network, events, timepoints, Vec::new())
    }

    /// Simulate composition with time-scheduled edge flow updates.
    pub fn simulate_with_flow_events<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut events: Vec<InletCompositionEvent<T>>,
        mut timepoints: Vec<T>,
        mut flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));
        flow_events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut active_flow_overrides: HashMap<usize, T> = HashMap::new();
        let mut event_cursor = 0usize;
        let mut flow_event_cursor = 0usize;
        let mut states = Vec::with_capacity(timepoints.len());

        for &time in &timepoints {
            while event_cursor < events.len() && events[event_cursor].time <= time {
                let event = &events[event_cursor];
                active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
                event_cursor += 1;
            }

            while flow_event_cursor < flow_events.len() && flow_events[flow_event_cursor].time <= time {
                let flow_event = &flow_events[flow_event_cursor];
                active_flow_overrides.insert(flow_event.edge_index, flow_event.flow_rate);
                flow_event_cursor += 1;
            }

            let mut effective_flow_rates: HashMap<usize, T> = HashMap::new();
            for (edge_idx, q) in &network.flow_rates {
                effective_flow_rates.insert(edge_idx.index(), *q);
            }
            for (edge_index, flow_rate) in &active_flow_overrides {
                effective_flow_rates.insert(*edge_index, *flow_rate);
            }

            let node_mixtures =
                Self::solve_node_mixtures(network, &active_inlet_mixtures, &effective_flow_rates)?;
            let edge_mixtures = Self::compute_edge_mixtures(network, &node_mixtures, &effective_flow_rates);

            states.push(CompositionState {
                time,
                node_mixtures,
                edge_mixtures,
                edge_flow_rates: effective_flow_rates,
            });
        }

        Ok(states)
    }

    /// Simulate composition with time-scheduled node pressure boundary updates.
    ///
    /// At each requested timepoint, pressure events active at that time are applied,
    /// the hydraulic network is re-solved, and composition mixing is computed from
    /// the resulting flow field.
    pub fn simulate_with_pressure_events<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        mut timepoints: Vec<T>,
        mut pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));
        pressure_events.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut working_network = network.clone();
        let solver = NetworkSolver::<T, F>::new();

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(timepoints.len());

        for &time in &timepoints {
            while composition_event_cursor < composition_events.len()
                && composition_events[composition_event_cursor].time <= time
            {
                let event = &composition_events[composition_event_cursor];
                active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
                composition_event_cursor += 1;
            }

            while pressure_event_cursor < pressure_events.len()
                && pressure_events[pressure_event_cursor].time <= time
            {
                let pressure_event = &pressure_events[pressure_event_cursor];
                working_network.set_pressure(
                    NodeIndex::new(pressure_event.node_index),
                    pressure_event.pressure,
                );
                pressure_event_cursor += 1;
            }

            let problem = NetworkProblem::new(working_network.clone());
            let solved_network = solver.solve_network(&problem)?;
            working_network = solved_network;

            let mut effective_flow_rates: HashMap<usize, T> = HashMap::new();
            for (edge_idx, q) in &working_network.flow_rates {
                effective_flow_rates.insert(edge_idx.index(), *q);
            }

            let node_mixtures = Self::solve_node_mixtures(
                &working_network,
                &active_inlet_mixtures,
                &effective_flow_rates,
            )?;
            let edge_mixtures =
                Self::compute_edge_mixtures(&working_network, &node_mixtures, &effective_flow_rates);

            states.push(CompositionState {
                time,
                node_mixtures,
                edge_mixtures,
                edge_flow_rates: effective_flow_rates,
            });
        }

        Ok(states)
    }

    fn times_close<T: RealField + Copy>(a: T, b: T, tolerance: T) -> bool {
        (a - b).abs() <= tolerance
    }

    fn sort_unique_timepoints<T: RealField + Copy>(
        mut timepoints: Vec<T>,
        tolerance: T,
    ) -> Vec<T> {
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut unique = Vec::with_capacity(timepoints.len());
        for time in timepoints {
            let should_push = unique
                .last()
                .map(|last| !Self::times_close(*last, time, tolerance))
                .unwrap_or(true);
            if should_push {
                unique.push(time);
            }
        }
        unique
    }

    fn solve_node_mixtures<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
    ) -> Result<HashMap<usize, MixtureComposition<T>>> {
        let mut node_mixtures = active_inlet_mixtures.clone();
        let max_iter = network.node_count().saturating_mul(4).max(8);
        let tolerance = T::from_f64(1e-9).unwrap_or_else(T::default_epsilon);

        for _ in 0..max_iter {
            let mut changed = false;

            for node_idx in network.graph.node_indices() {
                let node_id = node_idx.index();
                if active_inlet_mixtures.contains_key(&node_id) {
                    continue;
                }

                let mut incoming: Vec<(MixtureComposition<T>, T)> = Vec::new();

                for edge_ref in network.graph.edge_references() {
                    let src = edge_ref.source();
                    let dst = edge_ref.target();
                    let edge_idx = edge_ref.id();
                    let q = *flow_rates.get(&edge_idx.index()).unwrap_or(&T::zero());
                    let q_abs = q.abs();

                    if q_abs <= T::default_epsilon() {
                        continue;
                    }

                    if dst == node_idx && q > T::zero() {
                        if let Some(m) = node_mixtures.get(&src.index()) {
                            incoming.push((m.clone(), q_abs));
                        }
                    } else if src == node_idx && q < T::zero() {
                        if let Some(m) = node_mixtures.get(&dst.index()) {
                            incoming.push((m.clone(), q_abs));
                        }
                    }
                }

                if incoming.is_empty() {
                    continue;
                }

                let mixed = MixtureComposition::blend_weighted(&incoming);
                let should_update = match node_mixtures.get(&node_id) {
                    Some(current) => !current.approximately_equals(&mixed, tolerance),
                    None => true,
                };

                if should_update {
                    node_mixtures.insert(node_id, mixed);
                    changed = true;
                }
            }

            if !changed {
                return Ok(node_mixtures);
            }
        }

        Ok(node_mixtures)
    }

    fn compute_edge_mixtures<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
    ) -> HashMap<usize, MixtureComposition<T>> {
        let mut edge_mixtures = HashMap::new();

        for edge_ref in network.graph.edge_references() {
            let src = edge_ref.source();
            let dst = edge_ref.target();
            let edge_idx = edge_ref.id();
            let q = *flow_rates.get(&edge_idx.index()).unwrap_or(&T::zero());

            let upstream = if q >= T::zero() {
                src.index()
            } else {
                dst.index()
            };

            let composition = node_mixtures
                .get(&upstream)
                .cloned()
                .unwrap_or_else(MixtureComposition::empty);
            edge_mixtures.insert(edge_idx.index(), composition);
        }

        edge_mixtures
    }
}
