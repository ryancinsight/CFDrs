use super::events::{
    EdgeFlowEvent, InletCompositionEvent, InletHematocritEvent, PressureBoundaryEvent,
};
use super::state::{CompositionState, MixtureComposition};
use crate::domain::network::{
    Network, EDGE_PROPERTY_HEMATOCRIT, EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S,
    EDGE_PROPERTY_LOCAL_HEMATOCRIT, EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S,
};
use crate::solver::core::NetworkSolver;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

type CoupledBloodSnapshot<T, F> = (
    HashMap<usize, MixtureComposition<T>>,
    HashMap<usize, MixtureComposition<T>>,
    HashMap<usize, T>,
    Network<T, F>,
);

#[derive(Clone, Copy, Debug)]
struct IncidentEdge {
    edge_index: usize,
    source: usize,
    target: usize,
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
        Ok(Self::uniform_timepoints(
            self.duration,
            self.result_time_step,
        ))
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
        let tolerance = T::from_f64(1e-12).expect("Mathematical constant conversion compromised");
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
            .is_none_or(|last| (*last - duration).abs() > tolerance);

        if needs_terminal {
            points.push(duration);
        }

        points
    }
}

/// Simulator for transient inlet composition switching with instantaneous node mixing.
///
/// ## Algorithm
///
/// 1. **Event sorting**: All events (inlet composition changes, edge flow
///    overrides, and pressure boundary updates) are sorted by time and
///    merged with the uniform timepoint grid.
///
/// 2. **Time-stepping**: At each timepoint the simulator:
///    - Applies active `InletCompositionEvent`s to update inlet node mixtures.
///    - Applies active `EdgeFlowEvent`s to override edge flow rates.
///    - Applies active `PressureBoundaryEvent`s and re-solves the network
///      (pressure events trigger a full `NetworkSolver::solve_network` call).
///    - Propagates compositions downstream via flow-weighted instantaneous
///      mixing at each junction node.
///
/// 3. **Output**: A `Vec<CompositionState>` with one entry per timepoint,
///    recording the mixture composition at every node in the network.
pub struct TransientCompositionSimulator;

/// Configuration for segmented blood transport along each advecting edge.
#[derive(Debug, Clone)]
pub struct BloodEdgeTransportConfig<T: RealField + Copy> {
    /// Number of axial control volumes per transport-capable edge.
    pub segments_per_edge: usize,
    /// Maximum explicit Courant number used for substepping.
    pub max_courant_number: T,
}

impl<T: RealField + Copy + FromPrimitive> BloodEdgeTransportConfig<T> {
    /// Create a segmented transport configuration.
    #[must_use]
    pub fn new(segments_per_edge: usize, max_courant_number: T) -> Self {
        Self {
            segments_per_edge,
            max_courant_number,
        }
    }

    fn validate(&self) -> Result<()> {
        if self.segments_per_edge == 0 {
            return Err(Error::InvalidInput(
                "Segmented blood transport requires at least one segment per edge".to_string(),
            ));
        }
        if self.max_courant_number <= T::zero() {
            return Err(Error::InvalidInput(
                "Segmented blood transport requires a positive Courant limit".to_string(),
            ));
        }
        Ok(())
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for BloodEdgeTransportConfig<T> {
    fn default() -> Self {
        Self {
            segments_per_edge: 4,
            max_courant_number: T::from_f64(0.9)
                .expect("Mathematical constant conversion compromised"),
        }
    }
}

impl TransientCompositionSimulator {
    fn build_node_incidence_cache<T: RealField + Copy, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
    ) -> Vec<Vec<IncidentEdge>> {
        let node_count = network.node_count();
        let mut degree_counts = vec![0usize; node_count];

        for edge_ref in network.graph.edge_references() {
            let source = edge_ref.source().index();
            let target = edge_ref.target().index();
            degree_counts[source] += 1;
            if source != target {
                degree_counts[target] += 1;
            }
        }

        let mut incidence = degree_counts
            .into_iter()
            .map(Vec::with_capacity)
            .collect::<Vec<_>>();

        for edge_ref in network.graph.edge_references() {
            let edge = IncidentEdge {
                edge_index: edge_ref.id().index(),
                source: edge_ref.source().index(),
                target: edge_ref.target().index(),
            };
            incidence[edge.source].push(edge);
            if edge.source != edge.target {
                incidence[edge.target].push(edge);
            }
        }

        incidence
    }

    fn fill_edge_flow_rate_map<T: RealField + Copy, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge_flow_rates: &mut HashMap<usize, T>,
    ) {
        if edge_flow_rates.capacity() < network.edge_count() {
            edge_flow_rates.reserve(network.edge_count() - edge_flow_rates.capacity());
        }
        edge_flow_rates.clear();
        for (edge_index, &flow_rate) in network.flow_rates.iter().enumerate() {
            edge_flow_rates.insert(edge_index, flow_rate);
        }
    }

    fn fill_effective_flow_rates_from_overrides<T: RealField + Copy, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        active_flow_overrides: &HashMap<usize, T>,
        edge_flow_rates: &mut HashMap<usize, T>,
    ) {
        Self::fill_edge_flow_rate_map(network, edge_flow_rates);
        for (edge_index, flow_rate) in active_flow_overrides {
            edge_flow_rates.insert(*edge_index, *flow_rate);
        }
    }

    /// Simulate transported blood hematocrit using canonical RBC/plasma mixture keys.
    pub fn simulate_blood_hematocrit<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
        )
    }

    /// Simulate blood hematocrit with finite edge transport inventory.
    pub fn simulate_blood_hematocrit_with_edge_lag<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_hematocrit_with_edge_transport(network, events, timepoints)
    }

    /// Simulate blood hematocrit with a per-edge control-volume transport state.
    pub fn simulate_blood_hematocrit_with_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_mixture_with_flow_events_and_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            Vec::new(),
        )
    }

    /// Simulate blood hematocrit with segmented axial advection along each edge.
    pub fn simulate_blood_hematocrit_with_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        config.validate()?;
        Self::simulate_blood_mixture_with_flow_events_and_segmented_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            Vec::new(),
            config,
        )
    }

    /// Simulate transported blood hematocrit with time-scheduled edge-flow updates.
    pub fn simulate_blood_hematocrit_with_flow_events<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_with_flow_events(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            flow_events,
        )
    }

    /// Simulate blood hematocrit with edge-flow updates and finite edge transport.
    pub fn simulate_blood_hematocrit_with_flow_events_and_edge_lag<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_hematocrit_with_flow_events_and_edge_transport(
            network,
            events,
            timepoints,
            flow_events,
        )
    }

    /// Simulate blood hematocrit with edge-flow updates and per-edge transport state.
    pub fn simulate_blood_hematocrit_with_flow_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_mixture_with_flow_events_and_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            flow_events,
        )
    }

    /// Simulate blood hematocrit with flow events and segmented axial advection.
    pub fn simulate_blood_hematocrit_with_flow_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        flow_events: Vec<EdgeFlowEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        config.validate()?;
        Self::simulate_blood_mixture_with_flow_events_and_segmented_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            flow_events,
            config,
        )
    }

    /// Simulate transported blood hematocrit with pressure-boundary updates.
    pub fn simulate_blood_hematocrit_with_pressure_events<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_with_pressure_events(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
        )
    }

    /// Simulate blood hematocrit with pressure updates and per-edge transport state.
    pub fn simulate_blood_hematocrit_with_pressure_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_mixture_with_pressure_events_and_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
        )
    }

    /// Simulate blood hematocrit with pressure events and segmented axial advection.
    pub fn simulate_blood_hematocrit_with_pressure_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        config.validate()?;
        Self::simulate_blood_mixture_with_pressure_events_and_segmented_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
            config,
        )
    }

    /// Simulate transient blood hematocrit with pressure-event-driven hydraulic
    /// re-solves and hematocrit-dependent resistance coupling.
    pub fn simulate_blood_hematocrit_with_coupled_pressure_events<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_mixture_with_coupled_pressure_events(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
        )
    }

    /// Simulate transient blood hematocrit with coupled hydraulics and edge transport state.
    pub fn simulate_blood_hematocrit_with_coupled_pressure_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_mixture_with_coupled_pressure_events_and_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
        )
    }

    /// Simulate transient blood hematocrit with coupled hydraulics and segmented axial transport.
    pub fn simulate_blood_hematocrit_with_coupled_pressure_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timepoints: Vec<T>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        config.validate()?;
        Self::simulate_blood_mixture_with_coupled_pressure_events_and_segmented_edge_transport(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timepoints,
            pressure_events,
            config,
        )
    }

    /// Simulate transported blood hematocrit using MMFT-style timing controls.
    pub fn simulate_blood_hematocrit_with_time_config<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timing: SimulationTimeConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_with_time_config(
            network,
            Self::hematocrit_events_to_mixture_events(events)?,
            timing,
        )
    }

    /// Simulate blood hematocrit on MMFT-style timing controls with edge transport.
    pub fn simulate_blood_hematocrit_with_time_config_and_edge_lag<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timing: SimulationTimeConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        Self::simulate_blood_hematocrit_with_time_config_and_edge_transport(network, events, timing)
    }

    /// Simulate blood hematocrit on MMFT-style timing controls with edge transport state.
    pub fn simulate_blood_hematocrit_with_time_config_and_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timing: SimulationTimeConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        let result_timepoints = timing.result_timepoints()?;
        let mut internal_timepoints = result_timepoints.clone();
        internal_timepoints.extend(timing.calculation_timepoints()?);
        for event in &events {
            if event.time >= T::zero() && event.time <= timing.duration {
                internal_timepoints.push(event.time);
            }
        }
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);
        let states = Self::simulate_blood_hematocrit_with_edge_transport(
            network,
            events,
            merged_timepoints,
        )?;
        Self::sample_states_at_timepoints(states, result_timepoints, tolerance)
    }

    /// Simulate blood hematocrit on MMFT-style timing controls with segmented axial advection.
    pub fn simulate_blood_hematocrit_with_time_config_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        events: Vec<InletHematocritEvent<T>>,
        timing: SimulationTimeConfig<T>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        config.validate()?;
        let result_timepoints = timing.result_timepoints()?;
        let mut internal_timepoints = result_timepoints.clone();
        internal_timepoints.extend(timing.calculation_timepoints()?);
        for event in &events {
            if event.time >= T::zero() && event.time <= timing.duration {
                internal_timepoints.push(event.time);
            }
        }
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);
        let states = Self::simulate_blood_hematocrit_with_segmented_edge_transport(
            network,
            events,
            merged_timepoints,
            config,
        )?;
        Self::sample_states_at_timepoints(states, result_timepoints, tolerance)
    }

    fn simulate_blood_mixture_with_coupled_pressure_events<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
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

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut working_network = network.clone();
        let solver = NetworkSolver::<T, F>::new();
        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(timepoints.len());
        let tolerance = T::from_f64(1.0e-9).expect("Mathematical constant conversion compromised");
        let max_coupling_iters = working_network.node_count().clamp(2, 10);

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

            let mut current_flow_rates = HashMap::with_capacity(working_network.flow_rates.len());
            let mut current_node_mixtures = HashMap::with_capacity(working_network.node_count());
            let mut current_edge_mixtures =
                HashMap::with_capacity(working_network.flow_rates.len());
            let mut previous_flow_vector = working_network.flow_rates.clone();

            for _ in 0..max_coupling_iters {
                working_network = solver.solve_owned_network(working_network)?;

                current_flow_rates.clear();
                for (i, &q) in working_network.flow_rates.iter().enumerate() {
                    current_flow_rates.insert(i, q);
                }

                current_node_mixtures = Self::solve_node_mixtures(
                    &working_network,
                    &active_inlet_mixtures,
                    &current_flow_rates,
                )?;
                current_edge_mixtures = Self::compute_edge_mixtures(
                    &working_network,
                    &current_node_mixtures,
                    &current_flow_rates,
                );
                Self::backfill_blood_edge_mixtures_from_network(
                    &working_network,
                    &mut current_edge_mixtures,
                );

                let max_hct_change = Self::stamp_edge_hematocrit_from_mixtures(
                    &mut working_network,
                    &current_edge_mixtures,
                );
                let max_flow_change =
                    Self::max_flow_change(&previous_flow_vector, &working_network.flow_rates);

                if max_hct_change <= tolerance && max_flow_change <= tolerance {
                    break;
                }

                previous_flow_vector.clone_from(&working_network.flow_rates);
                working_network.update_resistances()?;
            }

            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures,
                edge_mixtures: current_edge_mixtures,
                edge_flow_rates: current_flow_rates,
            });
        }

        Ok(states)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_flow_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        flow_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(flow_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut active_flow_overrides: HashMap<usize, T> = HashMap::new();
        let mut event_cursor = 0usize;
        let mut flow_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_edge_mixtures = Self::initial_blood_edge_mixtures(network);
        let mut current_flow_rates = HashMap::with_capacity(network.edge_count());

        while event_cursor < composition_events.len()
            && composition_events[event_cursor].time <= previous_time
        {
            let event = &composition_events[event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            event_cursor += 1;
        }
        while flow_event_cursor < flow_events.len()
            && flow_events[flow_event_cursor].time <= previous_time
        {
            let event = &flow_events[flow_event_cursor];
            active_flow_overrides.insert(event.edge_index, event.flow_rate);
            flow_event_cursor += 1;
        }
        Self::fill_effective_flow_rates_from_overrides(
            network,
            &active_flow_overrides,
            &mut current_flow_rates,
        );
        let mut previous_flow_rates = current_flow_rates.clone();

        let mut current_node_mixtures = Self::solve_blood_node_mixtures_with_edge_transport(
            network,
            &active_inlet_mixtures,
            &previous_flow_rates,
            &previous_edge_mixtures,
        )?;
        let mut current_edge_snapshot = Self::compose_blood_edge_snapshot(
            network,
            &current_node_mixtures,
            &previous_flow_rates,
            &previous_edge_mixtures,
        );
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot.clone(),
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let transported_edge_mixtures = Self::advance_blood_edge_mixtures(
                network,
                &previous_edge_mixtures,
                &current_node_mixtures,
                &previous_flow_rates,
                dt,
            );

            while event_cursor < composition_events.len()
                && composition_events[event_cursor].time <= time
            {
                let event = &composition_events[event_cursor];
                active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
                event_cursor += 1;
            }
            while flow_event_cursor < flow_events.len()
                && flow_events[flow_event_cursor].time <= time
            {
                let event = &flow_events[flow_event_cursor];
                active_flow_overrides.insert(event.edge_index, event.flow_rate);
                flow_event_cursor += 1;
            }

            Self::fill_effective_flow_rates_from_overrides(
                network,
                &active_flow_overrides,
                &mut current_flow_rates,
            );
            current_node_mixtures = Self::solve_blood_node_mixtures_with_edge_transport(
                network,
                &active_inlet_mixtures,
                &current_flow_rates,
                &transported_edge_mixtures,
            )?;
            current_edge_snapshot = Self::compose_blood_edge_snapshot(
                network,
                &current_node_mixtures,
                &current_flow_rates,
                &transported_edge_mixtures,
            );
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot.clone(),
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates.clone_from(&current_flow_rates);
            previous_edge_mixtures = transported_edge_mixtures;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_flow_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut flow_events: Vec<EdgeFlowEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        flow_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(flow_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut active_flow_overrides: HashMap<usize, T> = HashMap::new();
        let mut event_cursor = 0usize;
        let mut flow_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_segment_state = Self::initial_segmented_blood_state(network, &config);
        let mut current_flow_rates = HashMap::with_capacity(network.edge_count());

        while event_cursor < composition_events.len()
            && composition_events[event_cursor].time <= previous_time
        {
            let event = &composition_events[event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            event_cursor += 1;
        }
        while flow_event_cursor < flow_events.len()
            && flow_events[flow_event_cursor].time <= previous_time
        {
            let event = &flow_events[flow_event_cursor];
            active_flow_overrides.insert(event.edge_index, event.flow_rate);
            flow_event_cursor += 1;
        }
        Self::fill_effective_flow_rates_from_overrides(
            network,
            &active_flow_overrides,
            &mut current_flow_rates,
        );
        let mut previous_flow_rates = current_flow_rates.clone();
        let mut current_node_mixtures =
            Self::solve_blood_node_mixtures_with_segmented_edge_transport(
                network,
                &active_inlet_mixtures,
                &previous_flow_rates,
                &previous_segment_state,
            )?;
        let mut current_edge_snapshot = Self::compose_segmented_blood_edge_snapshot(
            network,
            &current_node_mixtures,
            &previous_flow_rates,
            &previous_segment_state,
        );
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot.clone(),
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let current_edge_inlet_hematocrits =
                Self::solve_blood_edge_inlet_hematocrits_with_segmented_edge_transport(
                    network,
                    &current_node_mixtures,
                    &previous_flow_rates,
                )?;
            let transported_segment_state = Self::advance_blood_edge_segments(
                network,
                &previous_segment_state,
                &current_node_mixtures,
                &current_edge_inlet_hematocrits,
                &previous_flow_rates,
                dt,
                &config,
            );

            while event_cursor < composition_events.len()
                && composition_events[event_cursor].time <= time
            {
                let event = &composition_events[event_cursor];
                active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
                event_cursor += 1;
            }
            while flow_event_cursor < flow_events.len()
                && flow_events[flow_event_cursor].time <= time
            {
                let event = &flow_events[flow_event_cursor];
                active_flow_overrides.insert(event.edge_index, event.flow_rate);
                flow_event_cursor += 1;
            }

            Self::fill_effective_flow_rates_from_overrides(
                network,
                &active_flow_overrides,
                &mut current_flow_rates,
            );
            current_node_mixtures = Self::solve_blood_node_mixtures_with_segmented_edge_transport(
                network,
                &active_inlet_mixtures,
                &current_flow_rates,
                &transported_segment_state,
            )?;
            current_edge_snapshot = Self::compose_segmented_blood_edge_snapshot(
                network,
                &current_node_mixtures,
                &current_flow_rates,
                &transported_segment_state,
            );
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot.clone(),
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates.clone_from(&current_flow_rates);
            previous_segment_state = transported_segment_state;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_pressure_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(pressure_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let solver = NetworkSolver::<T, F>::new();
        let mut working_network = network.clone();
        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_edge_mixtures = Self::initial_blood_edge_mixtures(&working_network);

        while composition_event_cursor < composition_events.len()
            && composition_events[composition_event_cursor].time <= previous_time
        {
            let event = &composition_events[composition_event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            composition_event_cursor += 1;
        }
        while pressure_event_cursor < pressure_events.len()
            && pressure_events[pressure_event_cursor].time <= previous_time
        {
            let event = &pressure_events[pressure_event_cursor];
            working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
            pressure_event_cursor += 1;
        }
        working_network = solver.solve_owned_network(working_network)?;

        let mut previous_flow_rates = HashMap::with_capacity(working_network.edge_count());
        let mut current_flow_rates = HashMap::with_capacity(working_network.edge_count());
        Self::fill_edge_flow_rate_map(&working_network, &mut previous_flow_rates);
        let mut current_node_mixtures = Self::solve_blood_node_mixtures_with_edge_transport(
            &working_network,
            &active_inlet_mixtures,
            &previous_flow_rates,
            &previous_edge_mixtures,
        )?;
        let current_edge_snapshot = Self::compose_blood_edge_snapshot(
            &working_network,
            &current_node_mixtures,
            &previous_flow_rates,
            &previous_edge_mixtures,
        );
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot,
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let transported_edge_mixtures = Self::advance_blood_edge_mixtures(
                &working_network,
                &previous_edge_mixtures,
                &current_node_mixtures,
                &previous_flow_rates,
                dt,
            );

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
                let event = &pressure_events[pressure_event_cursor];
                working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
                pressure_event_cursor += 1;
            }

            working_network = solver.solve_owned_network(working_network)?;
            Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
            current_node_mixtures = Self::solve_blood_node_mixtures_with_edge_transport(
                &working_network,
                &active_inlet_mixtures,
                &current_flow_rates,
                &transported_edge_mixtures,
            )?;
            let current_edge_snapshot = Self::compose_blood_edge_snapshot(
                &working_network,
                &current_node_mixtures,
                &current_flow_rates,
                &transported_edge_mixtures,
            );
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot,
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates.clone_from(&current_flow_rates);
            previous_edge_mixtures = transported_edge_mixtures;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_pressure_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut pressure_events: Vec<PressureBoundaryEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(pressure_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let solver = NetworkSolver::<T, F>::new();
        let mut working_network = network.clone();
        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_segment_state =
            Self::initial_segmented_blood_state(&working_network, &config);

        while composition_event_cursor < composition_events.len()
            && composition_events[composition_event_cursor].time <= previous_time
        {
            let event = &composition_events[composition_event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            composition_event_cursor += 1;
        }
        while pressure_event_cursor < pressure_events.len()
            && pressure_events[pressure_event_cursor].time <= previous_time
        {
            let event = &pressure_events[pressure_event_cursor];
            working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
            pressure_event_cursor += 1;
        }
        working_network = solver.solve_owned_network(working_network)?;
        let mut previous_flow_rates = HashMap::with_capacity(working_network.edge_count());
        let mut current_flow_rates = HashMap::with_capacity(working_network.edge_count());
        Self::fill_edge_flow_rate_map(&working_network, &mut previous_flow_rates);
        let mut current_node_mixtures =
            Self::solve_blood_node_mixtures_with_segmented_edge_transport(
                &working_network,
                &active_inlet_mixtures,
                &previous_flow_rates,
                &previous_segment_state,
            )?;
        let mut current_edge_snapshot = Self::compose_segmented_blood_edge_snapshot(
            &working_network,
            &current_node_mixtures,
            &previous_flow_rates,
            &previous_segment_state,
        );
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot.clone(),
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let current_edge_inlet_hematocrits =
                Self::solve_blood_edge_inlet_hematocrits_with_segmented_edge_transport(
                    &working_network,
                    &current_node_mixtures,
                    &previous_flow_rates,
                )?;
            let transported_segment_state = Self::advance_blood_edge_segments(
                &working_network,
                &previous_segment_state,
                &current_node_mixtures,
                &current_edge_inlet_hematocrits,
                &previous_flow_rates,
                dt,
                &config,
            );

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
                let event = &pressure_events[pressure_event_cursor];
                working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
                pressure_event_cursor += 1;
            }

            working_network = solver.solve_owned_network(working_network)?;
            Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
            current_node_mixtures = Self::solve_blood_node_mixtures_with_segmented_edge_transport(
                &working_network,
                &active_inlet_mixtures,
                &current_flow_rates,
                &transported_segment_state,
            )?;
            current_edge_snapshot = Self::compose_segmented_blood_edge_snapshot(
                &working_network,
                &current_node_mixtures,
                &current_flow_rates,
                &transported_segment_state,
            );
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot.clone(),
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates.clone_from(&current_flow_rates);
            previous_segment_state = transported_segment_state;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_coupled_pressure_events_and_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(pressure_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let solver = NetworkSolver::<T, F>::new();
        let mut working_network = network.clone();
        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_edge_mixtures = Self::initial_blood_edge_mixtures(&working_network);
        let max_coupling_iters = working_network.node_count().clamp(2, 10);

        while composition_event_cursor < composition_events.len()
            && composition_events[composition_event_cursor].time <= previous_time
        {
            let event = &composition_events[composition_event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            composition_event_cursor += 1;
        }
        while pressure_event_cursor < pressure_events.len()
            && pressure_events[pressure_event_cursor].time <= previous_time
        {
            let event = &pressure_events[pressure_event_cursor];
            working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
            pressure_event_cursor += 1;
        }

        let (
            mut current_node_mixtures,
            mut current_edge_snapshot,
            mut previous_flow_rates,
            resolved_network,
        ) = Self::resolve_coupled_blood_snapshot(
            working_network,
            &active_inlet_mixtures,
            &previous_edge_mixtures,
            tolerance,
            max_coupling_iters,
            &solver,
        )?;
        working_network = resolved_network;
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot.clone(),
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let transported_edge_mixtures = Self::advance_blood_edge_mixtures(
                &working_network,
                &previous_edge_mixtures,
                &current_node_mixtures,
                &previous_flow_rates,
                dt,
            );

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
                let event = &pressure_events[pressure_event_cursor];
                working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
                pressure_event_cursor += 1;
            }

            let (resolved_nodes, resolved_edges, current_flow_rates, resolved_network) =
                Self::resolve_coupled_blood_snapshot(
                    working_network,
                    &active_inlet_mixtures,
                    &transported_edge_mixtures,
                    tolerance,
                    max_coupling_iters,
                    &solver,
                )?;
            working_network = resolved_network;
            current_node_mixtures = resolved_nodes;
            current_edge_snapshot = resolved_edges;
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot.clone(),
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates = current_flow_rates;
            previous_edge_mixtures = transported_edge_mixtures;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

    #[allow(clippy::too_many_lines)]
    fn simulate_blood_mixture_with_coupled_pressure_events_and_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        mut composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
        mut pressure_events: Vec<PressureBoundaryEvent<T>>,
        config: BloodEdgeTransportConfig<T>,
    ) -> Result<Vec<CompositionState<T>>> {
        if timepoints.is_empty() {
            return Err(Error::InvalidInput(
                "Transient composition simulation requires at least one timepoint".to_string(),
            ));
        }

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let requested_timepoints = Self::sorted_timepoints(timepoints);
        let mut internal_timepoints = requested_timepoints.clone();
        internal_timepoints.extend(composition_events.iter().map(|event| event.time));
        internal_timepoints.extend(pressure_events.iter().map(|event| event.time));
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(internal_timepoints, tolerance);

        let solver = NetworkSolver::<T, F>::new();
        let mut working_network = network.clone();
        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(merged_timepoints.len());
        let mut previous_time = merged_timepoints[0];
        let mut previous_segment_state =
            Self::initial_segmented_blood_state(&working_network, &config);
        let max_coupling_iters = working_network.node_count().clamp(2, 10);

        while composition_event_cursor < composition_events.len()
            && composition_events[composition_event_cursor].time <= previous_time
        {
            let event = &composition_events[composition_event_cursor];
            active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
            composition_event_cursor += 1;
        }
        while pressure_event_cursor < pressure_events.len()
            && pressure_events[pressure_event_cursor].time <= previous_time
        {
            let event = &pressure_events[pressure_event_cursor];
            working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
            pressure_event_cursor += 1;
        }

        let (
            mut current_node_mixtures,
            mut current_edge_snapshot,
            mut previous_flow_rates,
            resolved_network,
        ) = Self::resolve_coupled_segmented_blood_snapshot(
            working_network,
            &active_inlet_mixtures,
            &previous_segment_state,
            tolerance,
            max_coupling_iters,
            &solver,
        )?;
        working_network = resolved_network;
        states.push(CompositionState {
            time: previous_time,
            node_mixtures: current_node_mixtures.clone(),
            edge_mixtures: current_edge_snapshot.clone(),
            edge_flow_rates: previous_flow_rates.clone(),
        });

        for &time in merged_timepoints.iter().skip(1) {
            let dt = time - previous_time;
            let current_edge_inlet_hematocrits =
                Self::solve_blood_edge_inlet_hematocrits_with_segmented_edge_transport(
                    &working_network,
                    &current_node_mixtures,
                    &previous_flow_rates,
                )?;
            let transported_segment_state = Self::advance_blood_edge_segments(
                &working_network,
                &previous_segment_state,
                &current_node_mixtures,
                &current_edge_inlet_hematocrits,
                &previous_flow_rates,
                dt,
                &config,
            );

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
                let event = &pressure_events[pressure_event_cursor];
                working_network.set_pressure(NodeIndex::new(event.node_index), event.pressure);
                pressure_event_cursor += 1;
            }

            let (resolved_nodes, resolved_edges, current_flow_rates, resolved_network) =
                Self::resolve_coupled_segmented_blood_snapshot(
                    working_network,
                    &active_inlet_mixtures,
                    &transported_segment_state,
                    tolerance,
                    max_coupling_iters,
                    &solver,
                )?;
            working_network = resolved_network;
            current_node_mixtures = resolved_nodes;
            current_edge_snapshot = resolved_edges;
            states.push(CompositionState {
                time,
                node_mixtures: current_node_mixtures.clone(),
                edge_mixtures: current_edge_snapshot.clone(),
                edge_flow_rates: current_flow_rates.clone(),
            });

            previous_time = time;
            previous_flow_rates = current_flow_rates;
            previous_segment_state = transported_segment_state;
        }

        Self::sample_states_at_timepoints(states, requested_timepoints, tolerance)
    }

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

        let tolerance = T::from_f64(1e-12).expect("Mathematical constant conversion compromised");
        let merged_timepoints = Self::sort_unique_timepoints(all_timepoints, tolerance);
        let all_states = Self::simulate(network, events, merged_timepoints)?;
        Self::sample_states_at_timepoints(all_states, result_timepoints, tolerance).map_err(
            |source| Error::WithContext {
                context: "Failed to sample one or more configured result timepoints".to_string(),
                source: Box::new(source),
            },
        )
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

        events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        flow_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut active_flow_overrides: HashMap<usize, T> = HashMap::new();
        let mut event_cursor = 0usize;
        let mut flow_event_cursor = 0usize;
        let mut states = Vec::with_capacity(timepoints.len());
        let mut effective_flow_rates = HashMap::with_capacity(network.edge_count());

        for &time in &timepoints {
            while event_cursor < events.len() && events[event_cursor].time <= time {
                let event = &events[event_cursor];
                active_inlet_mixtures.insert(event.node_index, event.mixture.clone());
                event_cursor += 1;
            }

            while flow_event_cursor < flow_events.len()
                && flow_events[flow_event_cursor].time <= time
            {
                let flow_event = &flow_events[flow_event_cursor];
                active_flow_overrides.insert(flow_event.edge_index, flow_event.flow_rate);
                flow_event_cursor += 1;
            }

            Self::fill_effective_flow_rates_from_overrides(
                network,
                &active_flow_overrides,
                &mut effective_flow_rates,
            );

            let node_mixtures =
                Self::solve_node_mixtures(network, &active_inlet_mixtures, &effective_flow_rates)?;
            let edge_mixtures =
                Self::compute_edge_mixtures(network, &node_mixtures, &effective_flow_rates);

            states.push(CompositionState {
                time,
                node_mixtures,
                edge_mixtures,
                edge_flow_rates: effective_flow_rates.clone(),
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
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
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

        composition_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        pressure_events.sort_by(|a, b| {
            a.time
                .partial_cmp(&b.time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut working_network = network.clone();
        let solver = NetworkSolver::<T, F>::new();

        let mut active_inlet_mixtures: HashMap<usize, MixtureComposition<T>> = HashMap::new();
        let mut composition_event_cursor = 0usize;
        let mut pressure_event_cursor = 0usize;
        let mut states = Vec::with_capacity(timepoints.len());
        let mut effective_flow_rates = HashMap::with_capacity(network.edge_count());

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

            working_network = solver.solve_owned_network(working_network)?;
            Self::fill_edge_flow_rate_map(&working_network, &mut effective_flow_rates);

            let node_mixtures = Self::solve_node_mixtures(
                &working_network,
                &active_inlet_mixtures,
                &effective_flow_rates,
            )?;
            let edge_mixtures = Self::compute_edge_mixtures(
                &working_network,
                &node_mixtures,
                &effective_flow_rates,
            );

            states.push(CompositionState {
                time,
                node_mixtures,
                edge_mixtures,
                edge_flow_rates: effective_flow_rates.clone(),
            });
        }

        Ok(states)
    }

    fn times_close<T: RealField + Copy>(a: T, b: T, tolerance: T) -> bool {
        (a - b).abs() <= tolerance
    }

    fn hematocrit_events_to_mixture_events<T: RealField + Copy + FromPrimitive>(
        events: Vec<InletHematocritEvent<T>>,
    ) -> Result<Vec<InletCompositionEvent<T>>> {
        let mut converted = Vec::with_capacity(events.len());
        for event in events {
            if event.hematocrit < T::zero() || event.hematocrit > T::one() {
                return Err(Error::InvalidInput(
                    "Transient hematocrit events require values in [0, 1]".to_string(),
                ));
            }
            converted.push(InletCompositionEvent {
                time: event.time,
                node_index: event.node_index,
                mixture: MixtureComposition::from_blood_hematocrit(event.hematocrit),
            });
        }
        Ok(converted)
    }

    fn stamp_edge_hematocrit_from_mixtures<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &mut Network<T, F>,
        edge_mixtures: &HashMap<usize, MixtureComposition<T>>,
    ) -> T {
        let mut max_change = T::zero();
        for (edge_index, mixture) in edge_mixtures {
            if let Some(hematocrit) = mixture.hematocrit() {
                let edge_idx = petgraph::graph::EdgeIndex::new(*edge_index);
                if let Some(props) = network.properties.get_mut(&edge_idx) {
                    let current = props
                        .properties
                        .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                        .copied()
                        .unwrap_or(T::zero());
                    max_change = max_change.max((current - hematocrit).abs());
                    props.properties.insert(
                        EDGE_PROPERTY_LOCAL_HEMATOCRIT.to_string(),
                        hematocrit.max(T::zero()).min(T::one()),
                    );
                }
            }
        }
        max_change
    }

    fn stamp_edge_apparent_viscosity_from_segments<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &mut Network<T, F>,
        flow_rates: &HashMap<usize, T>,
        segment_state: &HashMap<usize, Vec<T>>,
    ) {
        for (edge_index, segments) in segment_state {
            let edge_idx = petgraph::graph::EdgeIndex::new(*edge_index);
            let Some(props) = network.properties.get_mut(&edge_idx) else {
                continue;
            };
            let plasma_viscosity = props
                .properties
                .get(EDGE_PROPERTY_PLASMA_VISCOSITY_PA_S)
                .copied();
            let flow_rate = flow_rates.get(edge_index).copied().unwrap_or_else(T::zero);
            let Some(segment_mu) = Self::segment_integrated_blood_apparent_viscosity(
                props.hydraulic_diameter,
                flow_rate,
                props.area,
                segments,
                plasma_viscosity,
            ) else {
                props
                    .properties
                    .remove(EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S);
                continue;
            };
            props.properties.insert(
                EDGE_PROPERTY_LOCAL_APPARENT_VISCOSITY_PA_S.to_string(),
                segment_mu,
            );
        }
    }

    fn max_flow_change<T: RealField + Copy>(previous: &[T], current: &[T]) -> T {
        let n = previous.len().max(current.len());
        let mut max_change = T::zero();
        for idx in 0..n {
            let a = previous.get(idx).copied().unwrap_or(T::zero());
            let b = current.get(idx).copied().unwrap_or(T::zero());
            max_change = max_change.max((a - b).abs());
        }
        max_change
    }

    fn backfill_blood_edge_mixtures_from_network<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        edge_mixtures: &mut HashMap<usize, MixtureComposition<T>>,
    ) {
        for (edge_idx, props) in &network.properties {
            let Some(hematocrit) = props
                .properties
                .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                .or_else(|| props.properties.get(EDGE_PROPERTY_HEMATOCRIT))
                .copied()
            else {
                continue;
            };

            let entry = edge_mixtures
                .entry(edge_idx.index())
                .or_insert_with(MixtureComposition::empty);
            if entry.hematocrit().is_none() {
                *entry = MixtureComposition::from_blood_hematocrit(hematocrit);
            }
        }
    }

    fn sort_unique_timepoints<T: RealField + Copy>(mut timepoints: Vec<T>, tolerance: T) -> Vec<T> {
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut unique = Vec::with_capacity(timepoints.len());
        for time in timepoints {
            let should_push = unique
                .last()
                .is_none_or(|last| !Self::times_close(*last, time, tolerance));
            if should_push {
                unique.push(time);
            }
        }
        unique
    }

    fn sorted_timepoints<T: RealField + Copy>(mut timepoints: Vec<T>) -> Vec<T> {
        timepoints.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        timepoints
    }

    fn sample_states_at_timepoints<T: RealField + Copy>(
        states: Vec<CompositionState<T>>,
        target_timepoints: Vec<T>,
        tolerance: T,
    ) -> Result<Vec<CompositionState<T>>> {
        let mut sampled: Vec<CompositionState<T>> = Vec::with_capacity(target_timepoints.len());
        let mut state_iter = states.into_iter();
        let mut current_state = state_iter.next();
        let mut last_matched_index: Option<usize> = None;

        for target_time in target_timepoints {
            if let Some(last_index) = last_matched_index {
                if Self::times_close(sampled[last_index].time, target_time, tolerance) {
                    sampled.push(sampled[last_index].clone());
                    continue;
                }
            }

            while let Some(state) = current_state.as_ref() {
                if state.time + tolerance < target_time {
                    current_state = state_iter.next();
                } else {
                    break;
                }
            }

            let Some(state) = current_state.take() else {
                return Err(Error::InvalidInput(
                    "Failed to sample one or more transient composition timepoints".to_string(),
                ));
            };
            if !Self::times_close(state.time, target_time, tolerance) {
                return Err(Error::InvalidInput(
                    "Failed to sample one or more transient composition timepoints".to_string(),
                ));
            }

            sampled.push(state);
            last_matched_index = Some(sampled.len() - 1);
            current_state = state_iter.next();
        }
        Ok(sampled)
    }

    fn initial_blood_edge_mixtures<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
    ) -> HashMap<usize, MixtureComposition<T>> {
        let mut edge_mixtures = HashMap::with_capacity(network.edge_count());
        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if Self::edge_supports_transport(network, edge_idx) {
                let hematocrit = network
                    .properties
                    .get(&edge_idx)
                    .and_then(|props| {
                        props
                            .properties
                            .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                            .or_else(|| props.properties.get(EDGE_PROPERTY_HEMATOCRIT))
                            .copied()
                    })
                    .unwrap_or_else(T::zero);
                edge_mixtures.insert(
                    edge_idx.index(),
                    MixtureComposition::from_blood_hematocrit(
                        hematocrit.max(T::zero()).min(T::one()),
                    ),
                );
            }
        }
        edge_mixtures
    }

    fn initial_segmented_blood_state<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        config: &BloodEdgeTransportConfig<T>,
    ) -> HashMap<usize, Vec<T>> {
        let mut segment_state = HashMap::with_capacity(network.edge_count());
        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if !Self::edge_supports_transport(network, edge_idx) {
                continue;
            }

            let hematocrit = network
                .properties
                .get(&edge_idx)
                .and_then(|props| {
                    props
                        .properties
                        .get(EDGE_PROPERTY_LOCAL_HEMATOCRIT)
                        .or_else(|| props.properties.get(EDGE_PROPERTY_HEMATOCRIT))
                        .copied()
                })
                .unwrap_or_else(T::zero);
            segment_state.insert(
                edge_idx.index(),
                vec![hematocrit.max(T::zero()).min(T::one()); config.segments_per_edge],
            );
        }
        segment_state
    }

    fn average_segment_hematocrit<T: RealField + Copy + FromPrimitive>(segments: &[T]) -> T {
        if segments.is_empty() {
            return T::zero();
        }
        let sum = segments
            .iter()
            .copied()
            .fold(T::zero(), |acc, value| acc + value);
        sum / T::from_usize(segments.len()).expect("usize conversion compromised")
    }

    fn outlet_segment_hematocrit<T: RealField + Copy + FromPrimitive>(
        segments: &[T],
        positive_direction: bool,
    ) -> T {
        if positive_direction {
            segments.last().copied().unwrap_or_else(T::zero)
        } else {
            segments.first().copied().unwrap_or_else(T::zero)
        }
    }

    fn segment_integrated_blood_apparent_viscosity<T: RealField + Copy + FromPrimitive>(
        hydraulic_diameter: Option<T>,
        flow_rate: T,
        area: T,
        segments: &[T],
        plasma_viscosity_override: Option<T>,
    ) -> Option<T> {
        let d_h = hydraulic_diameter?;
        if segments.is_empty() {
            return None;
        }
        let plasma_viscosity = plasma_viscosity_override.or_else(|| T::from_f64(1.05e-3))?;
        let mut mu_sum = T::zero();
        let mut count = 0usize;
        for &hematocrit in segments {
            let mu = crate::domain::network::blood_microchannel_apparent_viscosity(
                d_h,
                flow_rate,
                area,
                hematocrit,
                plasma_viscosity,
            )?;
            mu_sum += mu;
            count += 1;
        }
        if count == 0 {
            None
        } else {
            Some(mu_sum / T::from_usize(count).expect("usize conversion compromised"))
        }
    }

    fn resolve_coupled_blood_snapshot<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        mut working_network: Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        transported_edge_mixtures: &HashMap<usize, MixtureComposition<T>>,
        tolerance: T,
        max_coupling_iters: usize,
        solver: &NetworkSolver<T, F>,
    ) -> Result<CoupledBloodSnapshot<T, F>> {
        let mut current_node_mixtures = HashMap::new();
        let mut current_edge_snapshot = HashMap::new();
        let mut previous_flow_vector = working_network.flow_rates.clone();
        let mut current_flow_rates = HashMap::with_capacity(working_network.edge_count());
        Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);

        for _ in 0..max_coupling_iters {
            working_network = solver.solve_owned_network(working_network)?;
            Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
            current_node_mixtures = Self::solve_blood_node_mixtures_with_edge_transport(
                &working_network,
                active_inlet_mixtures,
                &current_flow_rates,
                transported_edge_mixtures,
            )?;
            current_edge_snapshot = Self::compose_blood_edge_snapshot(
                &working_network,
                &current_node_mixtures,
                &current_flow_rates,
                transported_edge_mixtures,
            );

            let max_hct_change = Self::stamp_edge_hematocrit_from_mixtures(
                &mut working_network,
                &current_edge_snapshot,
            );
            let max_flow_change =
                Self::max_flow_change(&previous_flow_vector, &working_network.flow_rates);
            if max_hct_change <= tolerance && max_flow_change <= tolerance {
                return Ok((
                    current_node_mixtures,
                    current_edge_snapshot,
                    current_flow_rates,
                    working_network,
                ));
            }

            previous_flow_vector.clone_from(&working_network.flow_rates);
            working_network.update_resistances()?;
        }

        Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
        Ok((
            current_node_mixtures,
            current_edge_snapshot,
            current_flow_rates,
            working_network,
        ))
    }

    fn resolve_coupled_segmented_blood_snapshot<
        T: RealField + Copy + FromPrimitive + ToPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        mut working_network: Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        segment_state: &HashMap<usize, Vec<T>>,
        tolerance: T,
        max_coupling_iters: usize,
        solver: &NetworkSolver<T, F>,
    ) -> Result<CoupledBloodSnapshot<T, F>> {
        let mut current_node_mixtures = HashMap::new();
        let mut current_edge_snapshot = HashMap::new();
        let mut previous_flow_vector = working_network.flow_rates.clone();
        let mut current_flow_rates = HashMap::with_capacity(working_network.edge_count());
        Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);

        for _ in 0..max_coupling_iters {
            working_network = solver.solve_owned_network(working_network)?;
            Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
            current_node_mixtures = Self::solve_blood_node_mixtures_with_segmented_edge_transport(
                &working_network,
                active_inlet_mixtures,
                &current_flow_rates,
                segment_state,
            )?;
            current_edge_snapshot = Self::compose_segmented_blood_edge_snapshot(
                &working_network,
                &current_node_mixtures,
                &current_flow_rates,
                segment_state,
            );
            Self::stamp_edge_apparent_viscosity_from_segments(
                &mut working_network,
                &current_flow_rates,
                segment_state,
            );

            let max_hct_change = Self::stamp_edge_hematocrit_from_mixtures(
                &mut working_network,
                &current_edge_snapshot,
            );
            let max_flow_change =
                Self::max_flow_change(&previous_flow_vector, &working_network.flow_rates);
            if max_hct_change <= tolerance && max_flow_change <= tolerance {
                return Ok((
                    current_node_mixtures,
                    current_edge_snapshot,
                    current_flow_rates,
                    working_network,
                ));
            }

            previous_flow_vector.clone_from(&working_network.flow_rates);
            working_network.update_resistances()?;
        }

        Self::fill_edge_flow_rate_map(&working_network, &mut current_flow_rates);
        Ok((
            current_node_mixtures,
            current_edge_snapshot,
            current_flow_rates,
            working_network,
        ))
    }

    fn edge_supports_transport<T: RealField + Copy, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge_idx: petgraph::graph::EdgeIndex,
    ) -> bool {
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        network
            .properties
            .get(&edge_idx)
            .is_some_and(|props| props.area > tolerance && props.length > tolerance)
    }

    fn advance_blood_edge_mixtures<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        previous_edge_mixtures: &HashMap<usize, MixtureComposition<T>>,
        previous_node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        previous_flow_rates: &HashMap<usize, T>,
        dt: T,
    ) -> HashMap<usize, MixtureComposition<T>> {
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let mut advanced = previous_edge_mixtures.clone();
        if dt <= tolerance {
            return advanced;
        }

        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if !Self::edge_supports_transport(network, edge_idx) {
                continue;
            }

            let flow_rate = previous_flow_rates
                .get(&edge_idx.index())
                .copied()
                .unwrap_or_else(T::zero);
            if num_traits::Float::abs(flow_rate) <= tolerance {
                continue;
            }

            let Some(props) = network.properties.get(&edge_idx) else {
                continue;
            };
            let residence_time = props.area * props.length / num_traits::Float::abs(flow_rate);
            if residence_time <= tolerance {
                continue;
            }

            let upstream_node = if flow_rate >= T::zero() {
                edge_ref.source().index()
            } else {
                edge_ref.target().index()
            };
            let inlet_hct = previous_node_mixtures
                .get(&upstream_node)
                .and_then(MixtureComposition::hematocrit)
                .unwrap_or_else(T::zero);
            let previous_hct = previous_edge_mixtures
                .get(&edge_idx.index())
                .and_then(MixtureComposition::hematocrit)
                .unwrap_or_else(T::zero);
            let decay = num_traits::Float::exp(-(dt / residence_time));
            let updated_hct = inlet_hct + (previous_hct - inlet_hct) * decay;
            advanced.insert(
                edge_idx.index(),
                MixtureComposition::from_blood_hematocrit(num_traits::Float::min(
                    num_traits::Float::max(updated_hct, T::zero()),
                    T::one(),
                )),
            );
        }

        advanced
    }

    fn advance_blood_edge_segments<
        T: RealField + Copy + FromPrimitive + num_traits::Float,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        previous_segment_state: &HashMap<usize, Vec<T>>,
        previous_node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        previous_edge_inlet_hematocrits: &HashMap<usize, T>,
        previous_flow_rates: &HashMap<usize, T>,
        dt: T,
        config: &BloodEdgeTransportConfig<T>,
    ) -> HashMap<usize, Vec<T>> {
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let mut advanced = previous_segment_state.clone();
        if dt <= tolerance {
            return advanced;
        }

        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if !Self::edge_supports_transport(network, edge_idx) {
                continue;
            }

            let flow_rate = previous_flow_rates
                .get(&edge_idx.index())
                .copied()
                .unwrap_or_else(T::zero);
            if num_traits::Float::abs(flow_rate) <= tolerance {
                continue;
            }

            let Some(props) = network.properties.get(&edge_idx) else {
                continue;
            };
            let Some(current_segments) = advanced.get_mut(&edge_idx.index()) else {
                continue;
            };
            if current_segments.is_empty() {
                continue;
            }

            let segment_count =
                T::from_usize(current_segments.len()).expect("usize conversion compromised");
            let segment_residence_time =
                props.area * props.length / (num_traits::Float::abs(flow_rate) * segment_count);
            if segment_residence_time <= tolerance {
                continue;
            }

            let upstream_node = if flow_rate >= T::zero() {
                edge_ref.source().index()
            } else {
                edge_ref.target().index()
            };
            let inlet_hct = previous_edge_inlet_hematocrits
                .get(&edge_idx.index())
                .copied()
                .or_else(|| {
                    previous_node_mixtures
                        .get(&upstream_node)
                        .and_then(MixtureComposition::hematocrit)
                })
                .unwrap_or_else(T::zero);

            let total_courant = dt / segment_residence_time;
            let steps = num_traits::Float::ceil(total_courant / config.max_courant_number)
                .to_usize()
                .unwrap_or(1)
                .max(1);
            let sub_dt = dt / T::from_usize(steps).expect("usize conversion compromised");
            let cfl = sub_dt / segment_residence_time;

            for _ in 0..steps {
                let previous = current_segments.clone();
                if flow_rate >= T::zero() {
                    current_segments[0] = previous[0] + cfl * (inlet_hct - previous[0]);
                    for i in 1..previous.len() {
                        current_segments[i] = previous[i] + cfl * (previous[i - 1] - previous[i]);
                    }
                } else {
                    let last = previous.len() - 1;
                    current_segments[last] = previous[last] + cfl * (inlet_hct - previous[last]);
                    for i in (0..last).rev() {
                        current_segments[i] = previous[i] + cfl * (previous[i + 1] - previous[i]);
                    }
                }

                for segment in current_segments.iter_mut() {
                    *segment = num_traits::Float::min(
                        num_traits::Float::max(*segment, T::zero()),
                        T::one(),
                    );
                }
            }
        }

        advanced
    }

    #[allow(clippy::too_many_lines)]
    fn solve_blood_edge_inlet_hematocrits_with_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive + ToPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
    ) -> Result<HashMap<usize, T>> {
        let tolerance = T::from_f64(1.0e-12).expect("Mathematical constant conversion compromised");
        let micron_scale = 1.0e6_f64;
        let tiny_fraction = 1.0e-15_f64;
        let to_f64 = |value: T| {
            value.to_f64().ok_or_else(|| {
                Error::InvalidConfiguration(
                    "Segmented blood transport requires finite scalar conversions".to_string(),
                )
            })
        };
        let mut edge_inlet_hematocrits = HashMap::with_capacity(network.edge_count());
        let node_incidence = Self::build_node_incidence_cache(network);
        let mut incoming: Vec<(usize, T)> = Vec::with_capacity(2);
        let mut outgoing: Vec<(usize, T)> = Vec::with_capacity(2);

        for node_idx in network.graph.node_indices() {
            let node_id = node_idx.index();
            let node_hct = node_mixtures
                .get(&node_id)
                .and_then(MixtureComposition::hematocrit)
                .unwrap_or_else(T::zero);

            incoming.clear();
            outgoing.clear();

            let incident_edges = &node_incidence[node_id];
            if incoming.capacity() < incident_edges.len() {
                incoming.reserve(incident_edges.len() - incoming.capacity());
            }
            if outgoing.capacity() < incident_edges.len() {
                outgoing.reserve(incident_edges.len() - outgoing.capacity());
            }

            for incident in incident_edges {
                let q = *flow_rates.get(&incident.edge_index).unwrap_or(&T::zero());
                let q_abs = q.abs();
                if q_abs <= T::default_epsilon() {
                    continue;
                }

                if incident.source == node_id {
                    if q > T::zero() {
                        outgoing.push((incident.edge_index, q_abs));
                    } else if q < T::zero() {
                        incoming.push((incident.edge_index, q_abs));
                    }
                } else if incident.target == node_id {
                    if q > T::zero() {
                        incoming.push((incident.edge_index, q_abs));
                    } else if q < T::zero() {
                        outgoing.push((incident.edge_index, q_abs));
                    }
                }
            }

            if outgoing.is_empty() {
                continue;
            }

            if incoming.len() == 1 && outgoing.len() == 2 {
                let parent_edge = incoming[0].0;
                let Some(parent_props) = network
                    .properties
                    .get(&petgraph::graph::EdgeIndex::new(parent_edge))
                else {
                    for &(edge_idx, _) in &outgoing {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };
                let Some(parent_diameter) = parent_props.hydraulic_diameter else {
                    for &(edge_idx, _) in &outgoing {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };

                let mut daughters = outgoing.clone();
                daughters.sort_by(|left, right| {
                    let left_d = network
                        .properties
                        .get(&petgraph::graph::EdgeIndex::new(left.0))
                        .and_then(|props| props.hydraulic_diameter)
                        .unwrap_or_else(T::zero);
                    let right_d = network
                        .properties
                        .get(&petgraph::graph::EdgeIndex::new(right.0))
                        .and_then(|props| props.hydraulic_diameter)
                        .unwrap_or_else(T::zero);
                    left_d
                        .partial_cmp(&right_d)
                        .unwrap_or(std::cmp::Ordering::Equal)
                        .then_with(|| left.0.cmp(&right.0))
                });

                let first = daughters[0];
                let second = daughters[1];
                let Some(first_props) = network
                    .properties
                    .get(&petgraph::graph::EdgeIndex::new(first.0))
                else {
                    for &(edge_idx, _) in &daughters {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };
                let Some(second_props) = network
                    .properties
                    .get(&petgraph::graph::EdgeIndex::new(second.0))
                else {
                    for &(edge_idx, _) in &daughters {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };
                let Some(first_diameter) = first_props.hydraulic_diameter else {
                    for &(edge_idx, _) in &daughters {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };
                let Some(second_diameter) = second_props.hydraulic_diameter else {
                    for &(edge_idx, _) in &daughters {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                };

                let q_first = first.1;
                let q_second = second.1;
                let total_q = q_first + q_second;
                if total_q <= tolerance {
                    for &(edge_idx, _) in &daughters {
                        edge_inlet_hematocrits.insert(edge_idx, node_hct);
                    }
                    continue;
                }

                let daughter_a =
                    crate::physics::cell_separation::plasma_skimming::pries_phase_separation(
                        to_f64(node_hct.max(T::zero()).min(T::one()))?,
                        to_f64((q_first / total_q).max(T::zero()).min(T::one()))?,
                        to_f64(first_diameter)? * micron_scale,
                        to_f64(second_diameter)? * micron_scale,
                        to_f64(parent_diameter)? * micron_scale,
                    )
                    .ok();

                if let Some(daughter_a) = daughter_a {
                    let q_second_fraction =
                        to_f64((q_second / total_q).max(T::zero()).min(T::one()))?;
                    let daughter_b_hct = if q_second_fraction > tiny_fraction {
                        T::from_f64(
                            ((1.0 - daughter_a.cell_fraction) * to_f64(node_hct)?
                                / q_second_fraction)
                                .clamp(0.0, 1.0),
                        )
                        .expect("Mathematical constant conversion compromised")
                    } else {
                        T::zero()
                    };
                    edge_inlet_hematocrits.insert(
                        first.0,
                        T::from_f64(daughter_a.daughter_hematocrit)
                            .expect("Mathematical constant conversion compromised"),
                    );
                    edge_inlet_hematocrits.insert(second.0, daughter_b_hct);
                    continue;
                }
            }

            for &(edge_idx, _) in &outgoing {
                edge_inlet_hematocrits.insert(edge_idx, node_hct);
            }
        }

        Ok(edge_inlet_hematocrits)
    }

    fn solve_blood_node_mixtures_with_edge_transport<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
        transported_edge_mixtures: &HashMap<usize, MixtureComposition<T>>,
    ) -> Result<HashMap<usize, MixtureComposition<T>>> {
        let mut node_mixtures = active_inlet_mixtures.clone();
        let max_iter = network.node_count().saturating_mul(4).max(8);
        let tolerance = T::from_f64(1e-9).expect("Mathematical constant conversion compromised");
        let node_incidence = Self::build_node_incidence_cache(network);
        let mut incoming: Vec<(MixtureComposition<T>, T)> = Vec::new();

        for _ in 0..max_iter {
            let mut changed = false;

            for node_idx in network.graph.node_indices() {
                let node_id = node_idx.index();
                if active_inlet_mixtures.contains_key(&node_id) {
                    continue;
                }

                incoming.clear();
                let incident_edges = &node_incidence[node_id];
                if incoming.capacity() < incident_edges.len() {
                    incoming.reserve(incident_edges.len() - incoming.capacity());
                }

                for incident in incident_edges {
                    let q = *flow_rates.get(&incident.edge_index).unwrap_or(&T::zero());
                    let q_abs = q.abs();

                    if q_abs <= T::default_epsilon() {
                        continue;
                    }

                    if incident.source == node_id {
                        if q < T::zero() {
                            if let Some(mixture) =
                                transported_edge_mixtures.get(&incident.edge_index)
                            {
                                incoming.push((mixture.clone(), q_abs));
                                continue;
                            }

                            if let Some(mixture) = node_mixtures.get(&incident.target) {
                                incoming.push((mixture.clone(), q_abs));
                            }
                        }
                    } else if incident.target == node_id && q > T::zero() {
                        if let Some(mixture) = transported_edge_mixtures.get(&incident.edge_index) {
                            incoming.push((mixture.clone(), q_abs));
                            continue;
                        }

                        if let Some(mixture) = node_mixtures.get(&incident.source) {
                            incoming.push((mixture.clone(), q_abs));
                        }
                    }
                }

                if incoming.is_empty() {
                    continue;
                }

                let mixed = MixtureComposition::blend_weighted_owned(&incoming);
                incoming.clear();

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

    fn solve_blood_node_mixtures_with_segmented_edge_transport<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
        segment_state: &HashMap<usize, Vec<T>>,
    ) -> Result<HashMap<usize, MixtureComposition<T>>> {
        let mut node_mixtures = active_inlet_mixtures.clone();
        let max_iter = network.node_count().saturating_mul(4).max(8);
        let tolerance = T::from_f64(1e-9).expect("Mathematical constant conversion compromised");
        let node_incidence = Self::build_node_incidence_cache(network);
        let mut incoming_owned: Vec<(MixtureComposition<T>, T)> = Vec::new();

        for _ in 0..max_iter {
            let mut changed = false;

            for node_idx in network.graph.node_indices() {
                let node_id = node_idx.index();
                if active_inlet_mixtures.contains_key(&node_id) {
                    continue;
                }

                incoming_owned.clear();
                let incident_edges = &node_incidence[node_id];
                if incoming_owned.capacity() < incident_edges.len() {
                    incoming_owned.reserve(incident_edges.len() - incoming_owned.capacity());
                }

                for incident in incident_edges {
                    let q = *flow_rates.get(&incident.edge_index).unwrap_or(&T::zero());
                    let q_abs = q.abs();

                    if q_abs <= T::default_epsilon() {
                        continue;
                    }

                    if incident.source == node_id {
                        if q < T::zero() {
                            if let Some(segments) = segment_state.get(&incident.edge_index) {
                                let edge_hct =
                                    Self::outlet_segment_hematocrit(segments, q >= T::zero());
                                incoming_owned.push((
                                    MixtureComposition::from_blood_hematocrit(edge_hct),
                                    q_abs,
                                ));
                                continue;
                            }

                            if let Some(mixture) = node_mixtures.get(&incident.target) {
                                incoming_owned.push((mixture.clone(), q_abs));
                            }
                        }
                    } else if incident.target == node_id && q > T::zero() {
                        if let Some(segments) = segment_state.get(&incident.edge_index) {
                            let edge_hct =
                                Self::outlet_segment_hematocrit(segments, q >= T::zero());
                            incoming_owned
                                .push((MixtureComposition::from_blood_hematocrit(edge_hct), q_abs));
                            continue;
                        }

                        if let Some(mixture) = node_mixtures.get(&incident.source) {
                            incoming_owned.push((mixture.clone(), q_abs));
                        }
                    }
                }

                if incoming_owned.is_empty() {
                    continue;
                }

                let mixed = MixtureComposition::blend_weighted_owned(&incoming_owned);
                incoming_owned.clear();

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

    fn compose_blood_edge_snapshot<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
        transported_edge_mixtures: &HashMap<usize, MixtureComposition<T>>,
    ) -> HashMap<usize, MixtureComposition<T>> {
        let mut edge_mixtures = HashMap::with_capacity(network.edge_count());
        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if let Some(mixture) = transported_edge_mixtures.get(&edge_idx.index()) {
                edge_mixtures.insert(edge_idx.index(), mixture.clone());
                continue;
            }

            let q = *flow_rates.get(&edge_idx.index()).unwrap_or(&T::zero());
            let upstream_node = if q >= T::zero() {
                edge_ref.source().index()
            } else {
                edge_ref.target().index()
            };
            let composition = node_mixtures
                .get(&upstream_node)
                .cloned()
                .unwrap_or_else(MixtureComposition::empty);
            edge_mixtures.insert(edge_idx.index(), composition);
        }
        edge_mixtures
    }

    fn compose_segmented_blood_edge_snapshot<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        node_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
        segment_state: &HashMap<usize, Vec<T>>,
    ) -> HashMap<usize, MixtureComposition<T>> {
        let mut edge_mixtures = HashMap::with_capacity(network.edge_count());
        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if let Some(segments) = segment_state.get(&edge_idx.index()) {
                edge_mixtures.insert(
                    edge_idx.index(),
                    MixtureComposition::from_blood_hematocrit(Self::average_segment_hematocrit(
                        segments,
                    )),
                );
                continue;
            }

            let q = *flow_rates.get(&edge_idx.index()).unwrap_or(&T::zero());
            let upstream_node = if q >= T::zero() {
                edge_ref.source().index()
            } else {
                edge_ref.target().index()
            };
            let composition = node_mixtures
                .get(&upstream_node)
                .cloned()
                .unwrap_or_else(MixtureComposition::empty);
            edge_mixtures.insert(edge_idx.index(), composition);
        }
        edge_mixtures
    }

    fn solve_node_mixtures<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        active_inlet_mixtures: &HashMap<usize, MixtureComposition<T>>,
        flow_rates: &HashMap<usize, T>,
    ) -> Result<HashMap<usize, MixtureComposition<T>>> {
        let mut node_mixtures = active_inlet_mixtures.clone();
        let max_iter = network.node_count().saturating_mul(4).max(8);
        let tolerance = T::from_f64(1e-9).expect("Mathematical constant conversion compromised");
        let node_incidence = Self::build_node_incidence_cache(network);
        let mut incoming: Vec<(MixtureComposition<T>, T)> = Vec::new();

        for _ in 0..max_iter {
            let mut changed = false;

            for node_idx in network.graph.node_indices() {
                let node_id = node_idx.index();
                if active_inlet_mixtures.contains_key(&node_id) {
                    continue;
                }

                incoming.clear();
                let incident_edges = &node_incidence[node_id];
                if incoming.capacity() < incident_edges.len() {
                    incoming.reserve(incident_edges.len() - incoming.capacity());
                }

                for incident in incident_edges {
                    let q = *flow_rates.get(&incident.edge_index).unwrap_or(&T::zero());
                    let q_abs = q.abs();

                    if q_abs <= T::default_epsilon() {
                        continue;
                    }

                    if incident.target == node_id && q > T::zero() {
                        if let Some(m) = node_mixtures.get(&incident.source) {
                            incoming.push((m.clone(), q_abs));
                        }
                    } else if incident.source == node_id && q < T::zero() {
                        if let Some(m) = node_mixtures.get(&incident.target) {
                            incoming.push((m.clone(), q_abs));
                        }
                    }
                }

                if incoming.is_empty() {
                    continue;
                }

                let mixed = MixtureComposition::blend_weighted_owned(&incoming);
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
        let mut edge_mixtures = HashMap::with_capacity(network.edge_count());

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
