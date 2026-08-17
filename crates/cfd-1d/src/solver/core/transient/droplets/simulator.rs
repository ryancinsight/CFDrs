use super::types::{
    ActiveDroplet, ChannelOccupancy, DropletBoundary, DropletBranch, DropletInjection,
    DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState, DropletTrackingState,
    SplitMode,
};
use crate::domain::network::{Network, NodeType};
use crate::solver::core::transient::composition::{
    CompositionState, EdgeFlowEvent, InletCompositionEvent, PressureBoundaryEvent,
    TransientCompositionSimulator,
};
use aequitas::systems::si::quantities::{Dimensionless, Time, Volume, VolumetricFlowRate};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement};
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Simulator for droplet state tracking and occupancy transitions.
///
/// ## Algorithm
///
/// 1. **Composition solve**: Compute steady-state or transient flow fields at
///    each timepoint using `TransientCompositionSimulator` (applies inlet
///    composition events, edge flow events, and pressure boundary events).
///
/// 2. **Injection activation**: At each timestep, activate scheduled
///    `DropletInjection` entries whose `inject_time` falls within the current
///    interval. Each injection places an `ActiveDroplet` at the head of the
///    specified edge.
///
/// 3. **Advection**: For each active droplet, advance its position along the
///    current edge by Δx = v · Δt where v = Q/A is the local mean velocity
///    from the composition-state flow field. When a droplet reaches the end
///    of an edge, it transitions to the downstream junction.
///
/// 4. **Junction routing**: At branching junctions, the `DropletSplitPolicy`
///    determines whether the droplet follows the highest-flow branch
///    (`FollowMaxFlow`), splits into child droplets (`SplitByFlowFraction`),
///    or is consumed.
///
/// 5. **Snapshot collection**: After each timestep, a `DropletTrackingState`
///    snapshot records every active droplet's position, velocity, and
///    channel occupancy for post-processing.
pub struct TransientDropletSimulator;

impl TransientDropletSimulator {
    /// Simulate droplet states with time-scheduled flow events.
    ///
    /// This convenience API first computes transient composition states using
    /// flow events and then runs droplet tracking with the default split policy.
    pub fn simulate_with_flow_events<
        T: CfdScalar + Copy + FloatElement,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<Time<T>>,
        flow_events: Vec<EdgeFlowEvent<T>>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        Self::simulate_with_flow_events_and_policy(
            network,
            injections,
            composition_events,
            timepoints,
            flow_events,
            DropletSplitPolicy::default(),
        )
    }

    /// Simulate droplet states with flow events and an explicit split policy.
    pub fn simulate_with_flow_events_and_policy<
        T: CfdScalar + Copy + FloatElement,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<Time<T>>,
        flow_events: Vec<EdgeFlowEvent<T>>,
        split_policy: DropletSplitPolicy<T>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        let composition_states = TransientCompositionSimulator::simulate_with_flow_events(
            network,
            composition_events,
            timepoints,
            flow_events,
        )?;

        Self::simulate_on_composition_with_policy(
            network,
            injections,
            composition_states,
            split_policy,
        )
    }

    /// Simulate droplet states with time-scheduled pressure boundary events.
    ///
    /// This convenience API first computes transient composition states using
    /// pressure events and then runs droplet tracking with the default split policy.
    pub fn simulate_with_pressure_events<T: CfdScalar, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<Time<T>>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        Self::simulate_with_pressure_events_and_policy(
            network,
            injections,
            composition_events,
            timepoints,
            pressure_events,
            DropletSplitPolicy::default(),
        )
    }

    /// Simulate droplet states with pressure events and an explicit split policy.
    pub fn simulate_with_pressure_events_and_policy<T: CfdScalar, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<Time<T>>,
        pressure_events: Vec<PressureBoundaryEvent<T>>,
        split_policy: DropletSplitPolicy<T>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        let composition_states = TransientCompositionSimulator::simulate_with_pressure_events(
            network,
            composition_events,
            timepoints,
            pressure_events,
        )?;

        Self::simulate_on_composition_with_policy(
            network,
            injections,
            composition_states,
            split_policy,
        )
    }

    /// Simulate droplet state transitions on top of transient composition states.
    pub fn simulate_on_composition<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_states: Vec<CompositionState<T>>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        Self::simulate_on_composition_with_policy(
            network,
            injections,
            composition_states,
            DropletSplitPolicy::default(),
        )
    }

    /// Simulate droplet states with an explicit split policy.
    #[allow(clippy::too_many_lines)]
    pub fn simulate_on_composition_with_policy<
        T: CfdScalar + Copy + FloatElement,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        mut composition_states: Vec<CompositionState<T>>,
        split_policy: DropletSplitPolicy<T>,
    ) -> Result<Vec<DropletTrackingState<T>>> {
        if composition_states.is_empty() {
            return Err(Error::InvalidInput(
                "Droplet simulation requires at least one composition state".to_string(),
            ));
        }

        composition_states.sort_by(|a, b| {
            a.time
                .into_base()
                .partial_cmp(&b.time.into_base())
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut injections_sorted = injections;
        injections_sorted.sort_by(|a, b| {
            a.injection_time
                .partial_cmp(&b.injection_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut active: HashMap<i32, ActiveDroplet<T>> = HashMap::new();
        let mut output = Vec::with_capacity(composition_states.len());

        let first_time = composition_states[0].time.into_base();
        let mut previous_time = first_time;
        let mut state_network = network.clone();

        for state in &composition_states {
            let state_time = state.time.into_base();
            let dt = (state_time - previous_time).max(T::ZERO);
            previous_time = state_time;
            for (edge_idx, flow_rate) in &state.edge_flow_rates {
                state_network.set_flow_rate(EdgeIndex::new(*edge_idx), *flow_rate);
            }

            for injection in &injections_sorted {
                let entry = active.entry(injection.droplet_id).or_insert(ActiveDroplet {
                    state: DropletState::Injection,
                    branches: Vec::new(),
                });

                if entry.state == DropletState::Injection
                    && state_time >= injection.injection_time.into_base()
                {
                    entry.state = DropletState::Network;
                    entry.branches = vec![DropletBranch {
                        channel_index: injection.channel_index,
                        center: Dimensionless::from_base(
                            injection
                                .relative_position
                                .into_base()
                                .max(T::ZERO)
                                .min(T::ONE),
                        ),
                        volume: injection.volume,
                    }];
                }
            }

            for injection in &injections_sorted {
                if let Some(droplet) = active.get_mut(&injection.droplet_id) {
                    if droplet.state == DropletState::Network {
                        Self::advance_droplet(&state_network, droplet, dt, &split_policy)?;
                        Self::merge_branches(&state_network, droplet)?;
                    }
                }
            }

            let mut snapshots = HashMap::new();
            for injection in &injections_sorted {
                let droplet = active.get(&injection.droplet_id).ok_or_else(|| {
                    Error::InvalidConfiguration("Missing droplet state".to_string())
                })?;

                let representative_position = droplet.branches.first().map(|b| DropletPosition {
                    channel_index: b.channel_index,
                    relative_position: b.center,
                });

                let local_mixture = representative_position
                    .as_ref()
                    .and_then(|pos| state.edge_mixtures.get(&pos.channel_index).cloned());

                let branch_count = droplet.branches.len();
                let mut occupancy_spans = Vec::with_capacity(branch_count);
                let mut boundaries = Vec::with_capacity(branch_count * 2);
                let mut total_volume = T::ZERO;

                if droplet.state == DropletState::Network {
                    for branch in &droplet.branches {
                        total_volume += branch.volume.into_base();
                        let (start, end) = Self::branch_interval(&state_network, branch)?;
                        occupancy_spans.push(ChannelOccupancy {
                            channel_index: branch.channel_index,
                            start: Dimensionless::from_base(start),
                            end: Dimensionless::from_base(end),
                        });
                        boundaries.push(DropletBoundary {
                            channel_index: branch.channel_index,
                            relative_position: Dimensionless::from_base(start),
                        });
                        boundaries.push(DropletBoundary {
                            channel_index: branch.channel_index,
                            relative_position: Dimensionless::from_base(end),
                        });
                    }
                }
                snapshots.insert(
                    injection.droplet_id,
                    DropletSnapshot {
                        droplet_id: injection.droplet_id,
                        state: droplet.state,
                        position: representative_position,
                        occupancy_spans,
                        boundaries,
                        total_volume: Volume::from_base(total_volume),
                        fluid_id: injection.fluid_id,
                        local_mixture,
                    },
                );
            }

            output.push(DropletTrackingState {
                time: Time::from_base(state_time),
                droplets: snapshots,
            });
        }

        Ok(output)
    }

    fn edge_area<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge: EdgeIndex,
    ) -> T {
        let one = T::ONE;
        network.properties.get(&edge).map_or(one, |p| {
            let area = p.area.into_base();
            if area > T::ZERO {
                area
            } else {
                one
            }
        })
    }

    fn edge_length<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge: EdgeIndex,
    ) -> T {
        let one = T::ONE;
        network.properties.get(&edge).map_or(one, |p| {
            let length = p.length.into_base();
            if length > T::ZERO {
                length
            } else {
                one
            }
        })
    }

    fn advance_droplet<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        droplet: &mut ActiveDroplet<T>,
        dt: T,
        split_policy: &DropletSplitPolicy<T>,
    ) -> Result<()> {
        if dt <= T::ZERO {
            return Ok(());
        }

        let mut new_branches = Vec::with_capacity(droplet.branches.len());
        let mut any_sink = false;
        let mut any_trapped = false;
        let old_branches = std::mem::take(&mut droplet.branches);

        for branch in old_branches {
            Self::advance_branch(
                network,
                branch,
                dt,
                &mut new_branches,
                &mut any_sink,
                &mut any_trapped,
                split_policy,
            )?;
        }

        droplet.branches = new_branches;
        if droplet.branches.is_empty() {
            droplet.state = if any_sink {
                DropletState::Sink
            } else {
                DropletState::Trapped
            };
        } else {
            droplet.state = DropletState::Network;
        }

        Ok(())
    }

    fn branch_interval<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        branch: &DropletBranch<T>,
    ) -> Result<(T, T)> {
        let (start_raw, end_raw) = Self::branch_interval_raw(network, branch)?;
        let start = start_raw.max(T::ZERO);
        let end = end_raw.min(T::ONE);
        Ok((start, end))
    }

    fn branch_interval_raw<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        branch: &DropletBranch<T>,
    ) -> Result<(T, T)> {
        let edge = EdgeIndex::new(branch.channel_index);
        let area = Self::edge_area(network, edge);
        let length = Self::edge_length(network, edge);
        let eps = <T as FloatElement>::from_f64(1e-12);
        if area <= eps || length <= eps {
            return Err(Error::InvalidConfiguration(
                "Edge area/length must be positive for finite-length droplet tracking".to_string(),
            ));
        }
        let half = T::ONE / (T::ONE + T::ONE);
        let frac_len = (branch.volume.into_base() / (area * length)).max(T::ZERO);
        let half_len = frac_len * half;
        let center = branch.center.into_base();
        let start = center - half_len;
        let end = center + half_len;
        Ok((start, end))
    }

    fn advance_branch<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        mut branch: DropletBranch<T>,
        dt: T,
        out_branches: &mut Vec<DropletBranch<T>>,
        any_sink: &mut bool,
        any_trapped: &mut bool,
        split_policy: &DropletSplitPolicy<T>,
    ) -> Result<()> {
        if dt <= T::ZERO {
            out_branches.push(branch);
            return Ok(());
        }

        let eps = <T as FloatElement>::from_f64(1e-12);
        let hops_remaining = network.edge_count().saturating_mul(4).max(8);

        if let Some(_hop) = (0..hops_remaining).next() {
            let edge_idx = EdgeIndex::new(branch.channel_index);
            let q = network
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .map_or(T::ZERO, VolumetricFlowRate::into_base);
            if <T as NumericElement>::abs(q) <= eps {
                out_branches.push(branch);
                return Ok(());
            }

            let area = Self::edge_area(network, edge_idx);
            let length = Self::edge_length(network, edge_idx);
            if area <= eps || length <= eps {
                return Err(Error::InvalidConfiguration(
                    "Edge area/length must be positive for droplet advection".to_string(),
                ));
            }

            let dp = dt * (q / area) / length;
            let new_center = branch.center.into_base() + dp;
            let moved_branch = DropletBranch {
                channel_index: branch.channel_index,
                center: Dimensionless::from_base(new_center),
                volume: branch.volume,
            };
            let (start2, end) = Self::branch_interval_raw(network, &moved_branch)?;

            let crosses_downstream =
                (q >= T::ZERO && end > T::ONE - eps) || (q < T::ZERO && start2 < T::ZERO + eps);

            if !crosses_downstream {
                branch.center = Dimensionless::from_base(new_center.max(T::ZERO).min(T::ONE));
                out_branches.push(branch);
                return Ok(());
            }

            let (src, dst) = network
                .graph
                .edge_endpoints(edge_idx)
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge endpoints".to_string()))?;
            let boundary_node = if q >= T::ZERO { dst } else { src };

            if Self::is_sink_node(network, boundary_node)? {
                *any_sink = true;
                return Ok(());
            }

            let outgoing = Self::select_outgoing_edges(network, boundary_node, edge_idx);
            if outgoing.is_empty() {
                *any_trapped = true;
                return Ok(());
            }

            let selected = Self::select_split_targets(branch.volume, outgoing, split_policy, eps);
            if selected.is_empty() {
                *any_trapped = true;
                return Ok(());
            }

            let total_flow = selected
                .iter()
                .map(|(_, _, flow)| *flow)
                .fold(T::ZERO, |acc, v| acc + v);
            if total_flow <= eps {
                *any_trapped = true;
                return Ok(());
            }

            // Flow-weighted split with exact volume conservation.
            for (next_edge, start_position, flow_mag) in selected {
                let fraction = flow_mag / total_flow;
                let child_volume = Volume::from_base(branch.volume.into_base() * fraction);
                if child_volume.into_base() <= eps {
                    continue;
                }

                let next_area = Self::edge_area(network, next_edge);
                let next_length = Self::edge_length(network, next_edge);
                if next_area <= eps || next_length <= eps {
                    continue;
                }

                let frac_len = child_volume.into_base() / (next_area * next_length);
                let half = T::ONE / (T::ONE + T::ONE);
                let half_len = frac_len * half;
                let center = if start_position <= T::ZERO + eps {
                    half_len
                } else {
                    T::ONE - half_len
                }
                .max(T::ZERO)
                .min(T::ONE);

                out_branches.push(DropletBranch {
                    channel_index: next_edge.index(),
                    center: Dimensionless::from_base(center),
                    volume: child_volume,
                });
            }
            return Ok(());
        }

        *any_trapped = true;
        Ok(())
    }

    fn select_split_targets<T: CfdScalar + Copy + FloatElement>(
        branch_volume: Volume<T>,
        mut outgoing: Vec<(EdgeIndex, T, T)>,
        split_policy: &DropletSplitPolicy<T>,
        eps: T,
    ) -> Vec<(EdgeIndex, T, T)> {
        if outgoing.is_empty() {
            return Vec::new();
        }

        outgoing.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        let max_branches = split_policy.max_split_branches.max(1);
        if outgoing.len() > max_branches {
            outgoing.truncate(max_branches);
        }

        let dominant_only = vec![outgoing[0]];
        if outgoing.len() == 1 {
            return dominant_only;
        }

        match split_policy.mode {
            SplitMode::NeverSplit => dominant_only,
            SplitMode::AlwaysSplit => {
                let total = outgoing.iter().map(|x| x.2).fold(T::ZERO, |a, b| a + b);
                if total <= eps {
                    return dominant_only;
                }
                let mut min_child = branch_volume.into_base();
                for (_, _, flow_mag) in &outgoing {
                    let child = branch_volume.into_base() * (*flow_mag / total);
                    if child < min_child {
                        min_child = child;
                    }
                }
                if min_child < split_policy.min_child_volume.into_base() {
                    dominant_only
                } else {
                    outgoing
                }
            }
            SplitMode::AutoFlowWeighted => {
                let total = outgoing.iter().map(|x| x.2).fold(T::ZERO, |a, b| a + b);
                if total <= eps {
                    return dominant_only;
                }

                let dominant = outgoing[0].2;
                let secondary = total - dominant;
                let secondary_frac = secondary / total;

                if secondary_frac < split_policy.min_secondary_flow_fraction.into_base() {
                    return dominant_only;
                }

                let mut min_child = branch_volume.into_base();
                for (_, _, flow_mag) in &outgoing {
                    let child = branch_volume.into_base() * (*flow_mag / total);
                    if child < min_child {
                        min_child = child;
                    }
                }
                if min_child < split_policy.min_child_volume.into_base() {
                    dominant_only
                } else {
                    outgoing
                }
            }
        }
    }

    fn merge_branches<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        droplet: &mut ActiveDroplet<T>,
    ) -> Result<()> {
        if droplet.branches.len() <= 1 {
            return Ok(());
        }

        // Merge overlapping segments on same channel (volume-conserving).
        let mut groups: HashMap<usize, Vec<DropletBranch<T>>> = HashMap::new();
        for b in droplet.branches.drain(..) {
            groups.entry(b.channel_index).or_default().push(b);
        }

        let mut merged = Vec::new();
        for (_channel, mut branches) in groups {
            branches.sort_by(|a, b| {
                a.center
                    .partial_cmp(&b.center)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut current = branches[0].clone();
            for next in branches.into_iter().skip(1) {
                let (c_start, c_end) = Self::branch_interval(network, &current)?;
                let (n_start, n_end) = Self::branch_interval(network, &next)?;
                if n_start <= c_end {
                    // Overlap -> merge by volume-weighted center.
                    let v_total = current.volume.into_base() + next.volume.into_base();
                    let center = if v_total > T::ZERO {
                        (current.center.into_base() * current.volume.into_base()
                            + next.center.into_base() * next.volume.into_base())
                            / v_total
                    } else {
                        current.center.into_base()
                    };
                    current = DropletBranch {
                        channel_index: current.channel_index,
                        center: Dimensionless::from_base(
                            center.max(c_start.min(n_start)).min(c_end.max(n_end)),
                        ),
                        volume: Volume::from_base(v_total),
                    };
                } else {
                    merged.push(current);
                    current = next;
                }
            }
            merged.push(current);
        }
        droplet.branches = merged;
        Ok(())
    }

    fn is_sink_node<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        node: NodeIndex,
    ) -> Result<bool> {
        let node_weight = network
            .graph
            .node_weight(node)
            .ok_or_else(|| Error::InvalidConfiguration("Node not found".to_string()))?;
        if matches!(node_weight.node_type, NodeType::Outlet) {
            return Ok(true);
        }

        // Treat terminal dead-end nodes as sinks as well.
        let degree = network.graph.neighbors_undirected(node).count();
        Ok(degree <= 1)
    }

    fn select_outgoing_edges<T: CfdScalar + Copy + FloatElement, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        node: NodeIndex,
        previous_edge: EdgeIndex,
    ) -> Vec<(EdgeIndex, T, T)> {
        let eps = <T as FloatElement>::from_f64(1e-12);
        let mut candidates = Vec::new();

        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if edge_idx == previous_edge {
                continue;
            }

            let src = edge_ref.source();
            let dst = edge_ref.target();
            let q = network
                .flow_rates
                .get(edge_idx.index())
                .copied()
                .map_or(T::ZERO, VolumetricFlowRate::into_base);

            if <T as NumericElement>::abs(q) <= eps {
                continue;
            }

            let candidate = if src == node && q > T::ZERO {
                Some((edge_idx, T::ZERO, <T as NumericElement>::abs(q)))
            } else if dst == node && q < T::ZERO {
                Some((edge_idx, T::ONE, <T as NumericElement>::abs(q)))
            } else {
                None
            };

            if let Some((idx, start, mag)) = candidate {
                candidates.push((idx, start, mag));
            }
        }

        candidates
    }
}
