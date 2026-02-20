use super::types::{
    ActiveDroplet, ChannelOccupancy, DropletBoundary, DropletBranch, DropletInjection,
    DropletPosition, DropletSnapshot, DropletSplitPolicy, DropletState, DropletTrackingState,
    SplitMode,
};
use crate::solver::transient::composition::{
    CompositionState, EdgeFlowEvent, InletCompositionEvent, PressureBoundaryEvent,
    TransientCompositionSimulator,
};
use crate::network::{Network, NodeType};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;

/// Simulator for droplet state tracking and occupancy transitions.
pub struct TransientDropletSimulator;

impl TransientDropletSimulator {
    /// Simulate droplet states with time-scheduled flow events.
    ///
    /// This convenience API first computes transient composition states using
    /// flow events and then runs droplet tracking with the default split policy.
    pub fn simulate_with_flow_events<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
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
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
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
    pub fn simulate_with_pressure_events<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
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
    pub fn simulate_with_pressure_events_and_policy<
        T: RealField + Copy + FromPrimitive,
        F: FluidTrait<T> + Clone,
    >(
        network: &Network<T, F>,
        injections: Vec<DropletInjection<T>>,
        composition_events: Vec<InletCompositionEvent<T>>,
        timepoints: Vec<T>,
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
    pub fn simulate_on_composition<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
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
    pub fn simulate_on_composition_with_policy<
        T: RealField + Copy + FromPrimitive,
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

        composition_states.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal));

        let mut injections_sorted = injections;
        injections_sorted.sort_by(|a, b| {
            a.injection_time
                .partial_cmp(&b.injection_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut active: HashMap<i32, ActiveDroplet<T>> = HashMap::new();
        let mut output = Vec::with_capacity(composition_states.len());

        let first_time = composition_states[0].time;
        let mut previous_time = first_time;

        for state in &composition_states {
            let dt = (state.time - previous_time).max(T::zero());
            previous_time = state.time;

            let mut state_network = network.clone();
            for (edge_idx, flow_rate) in &state.edge_flow_rates {
                state_network.set_flow_rate(EdgeIndex::new(*edge_idx), *flow_rate);
            }

            for injection in &injections_sorted {
                let entry = active.entry(injection.droplet_id).or_insert(ActiveDroplet {
                    state: DropletState::Injection,
                    branches: Vec::new(),
                });

                if entry.state == DropletState::Injection && state.time >= injection.injection_time {
                    entry.state = DropletState::Network;
                    entry.branches = vec![DropletBranch {
                        channel_index: injection.channel_index,
                        center: injection.relative_position.max(T::zero()).min(T::one()),
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
                let droplet = active
                    .get(&injection.droplet_id)
                    .ok_or_else(|| Error::InvalidConfiguration("Missing droplet state".to_string()))?;

                let representative_position = droplet.branches.first().map(|b| DropletPosition {
                    channel_index: b.channel_index,
                    relative_position: b.center,
                });

                let local_mixture = representative_position
                    .as_ref()
                    .and_then(|pos| state.edge_mixtures.get(&pos.channel_index).cloned());

                let mut occupied_channels = Vec::new();
                let mut occupancy_spans = Vec::new();
                let mut boundaries = Vec::new();
                let mut total_volume = T::zero();

                if droplet.state == DropletState::Network {
                    for branch in &droplet.branches {
                        if !occupied_channels.contains(&branch.channel_index) {
                            occupied_channels.push(branch.channel_index);
                        }
                        total_volume += branch.volume;
                        let (start, end) = Self::branch_interval(&state_network, branch)?;
                        occupancy_spans.push(ChannelOccupancy {
                            channel_index: branch.channel_index,
                            start,
                            end,
                        });
                        boundaries.push(DropletBoundary {
                            channel_index: branch.channel_index,
                            relative_position: start,
                        });
                        boundaries.push(DropletBoundary {
                            channel_index: branch.channel_index,
                            relative_position: end,
                        });
                    }
                }

                snapshots.insert(
                    injection.droplet_id,
                    DropletSnapshot {
                        droplet_id: injection.droplet_id,
                        state: droplet.state,
                        position: representative_position,
                        occupied_channels,
                        occupancy_spans,
                        boundaries,
                        total_volume,
                        fluid_id: injection.fluid_id,
                        local_mixture,
                    },
                );
            }

            output.push(DropletTrackingState {
                time: state.time,
                droplets: snapshots,
            });
        }

        Ok(output)
    }

    fn edge_area<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge: EdgeIndex,
    ) -> T {
        let one = T::one();
        network
            .properties
            .get(&edge)
            .map(|p| if p.area > T::zero() { p.area } else { one })
            .unwrap_or(one)
    }

    fn edge_length<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        edge: EdgeIndex,
    ) -> T {
        let one = T::one();
        network
            .properties
            .get(&edge)
            .map(|p| if p.length > T::zero() { p.length } else { one })
            .unwrap_or(one)
    }

    fn advance_droplet<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        droplet: &mut ActiveDroplet<T>,
        dt: T,
        split_policy: &DropletSplitPolicy<T>,
    ) -> Result<()> {
        if dt <= T::zero() {
            return Ok(());
        }

        let mut new_branches = Vec::new();
        let mut any_sink = false;
        let mut any_trapped = false;

        for branch in droplet.branches.clone() {
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
            } else if any_trapped {
                DropletState::Trapped
            } else {
                DropletState::Trapped
            };
        } else {
            droplet.state = DropletState::Network;
        }

        Ok(())
    }

    fn branch_interval<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        branch: &DropletBranch<T>,
    ) -> Result<(T, T)> {
        let (start_raw, end_raw) = Self::branch_interval_raw(network, branch)?;
        let start = start_raw.max(T::zero());
        let end = end_raw.min(T::one());
        Ok((start, end))
    }

    fn branch_interval_raw<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        branch: &DropletBranch<T>,
    ) -> Result<(T, T)> {
        let edge = EdgeIndex::new(branch.channel_index);
        let area = Self::edge_area(network, edge);
        let length = Self::edge_length(network, edge);
        let eps = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);
        if area <= eps || length <= eps {
            return Err(Error::InvalidConfiguration(
                "Edge area/length must be positive for finite-length droplet tracking".to_string(),
            ));
        }
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        let frac_len = (branch.volume / (area * length)).max(T::zero());
        let half_len = frac_len * half;
        let start = branch.center - half_len;
        let end = branch.center + half_len;
        Ok((start, end))
    }

    fn advance_branch<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        mut branch: DropletBranch<T>,
        dt: T,
        out_branches: &mut Vec<DropletBranch<T>>,
        any_sink: &mut bool,
        any_trapped: &mut bool,
        split_policy: &DropletSplitPolicy<T>,
    ) -> Result<()> {
        if dt <= T::zero() {
            out_branches.push(branch);
            return Ok(());
        }

        let eps = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);
        let hops_remaining = network.edge_count().saturating_mul(4).max(8);

        for _hop in 0..hops_remaining {

            let edge_idx = EdgeIndex::new(branch.channel_index);
            let q = *network.flow_rates.get(&edge_idx).unwrap_or(&T::zero());
            if q.abs() <= eps {
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
            let new_center = branch.center + dp;
            let (start2, end) = Self::branch_interval_raw(
                network,
                &DropletBranch {
                    center: new_center,
                    ..branch.clone()
                },
            )?;

            let crosses_downstream = (q >= T::zero() && end > T::one() - eps)
                || (q < T::zero() && start2 < T::zero() + eps);

            if !crosses_downstream {
                branch.center = new_center.max(T::zero()).min(T::one());
                out_branches.push(branch);
                return Ok(());
            }

            let (src, dst) = network
                .graph
                .edge_endpoints(edge_idx)
                .ok_or_else(|| Error::InvalidConfiguration("Missing edge endpoints".to_string()))?;
            let boundary_node = if q >= T::zero() { dst } else { src };

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
                .fold(T::zero(), |acc, v| acc + v);
            if total_flow <= eps {
                *any_trapped = true;
                return Ok(());
            }

            // Flow-weighted split with exact volume conservation.
            for (next_edge, start_position, flow_mag) in selected {
                let fraction = flow_mag / total_flow;
                let child_volume = branch.volume * fraction;
                if child_volume <= eps {
                    continue;
                }

                let next_area = Self::edge_area(network, next_edge);
                let next_length = Self::edge_length(network, next_edge);
                if next_area <= eps || next_length <= eps {
                    continue;
                }

                let frac_len = child_volume / (next_area * next_length);
                let half = T::from_f64(0.5).unwrap_or_else(T::one);
                let half_len = frac_len * half;
                let center = if start_position <= T::zero() + eps {
                    half_len
                } else {
                    T::one() - half_len
                }
                .max(T::zero())
                .min(T::one());

                out_branches.push(DropletBranch {
                    channel_index: next_edge.index(),
                    center,
                    volume: child_volume,
                });
            }
            return Ok(());
        }

        *any_trapped = true;
        Ok(())
    }

    fn select_split_targets<T: RealField + Copy + FromPrimitive>(
        branch_volume: T,
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
                let total = outgoing.iter().map(|x| x.2).fold(T::zero(), |a, b| a + b);
                if total <= eps {
                    return dominant_only;
                }
                let mut min_child = branch_volume;
                for (_, _, flow_mag) in &outgoing {
                    let child = branch_volume * (*flow_mag / total);
                    if child < min_child {
                        min_child = child;
                    }
                }
                if min_child < split_policy.min_child_volume {
                    dominant_only
                } else {
                    outgoing
                }
            }
            SplitMode::AutoFlowWeighted => {
                let total = outgoing.iter().map(|x| x.2).fold(T::zero(), |a, b| a + b);
                if total <= eps {
                    return dominant_only;
                }

                let dominant = outgoing[0].2;
                let secondary = total - dominant;
                let secondary_frac = secondary / total;

                if secondary_frac < split_policy.min_secondary_flow_fraction {
                    return dominant_only;
                }

                let mut min_child = branch_volume;
                for (_, _, flow_mag) in &outgoing {
                    let child = branch_volume * (*flow_mag / total);
                    if child < min_child {
                        min_child = child;
                    }
                }
                if min_child < split_policy.min_child_volume {
                    dominant_only
                } else {
                    outgoing
                }
            }
        }
    }

    fn merge_branches<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
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
            branches.sort_by(|a, b| a.center.partial_cmp(&b.center).unwrap_or(std::cmp::Ordering::Equal));
            let mut current = branches[0].clone();
            for next in branches.into_iter().skip(1) {
                let (c_start, c_end) = Self::branch_interval(network, &current)?;
                let (n_start, n_end) = Self::branch_interval(network, &next)?;
                if n_start <= c_end {
                    // Overlap -> merge by volume-weighted center.
                    let v_total = current.volume + next.volume;
                    let center = if v_total > T::zero() {
                        (current.center * current.volume + next.center * next.volume) / v_total
                    } else {
                        current.center
                    };
                    current = DropletBranch {
                        channel_index: current.channel_index,
                        center: center.max(c_start.min(n_start)).min(c_end.max(n_end)),
                        volume: v_total,
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

    fn is_sink_node<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
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

    fn select_outgoing_edges<T: RealField + Copy + FromPrimitive, F: FluidTrait<T> + Clone>(
        network: &Network<T, F>,
        node: NodeIndex,
        previous_edge: EdgeIndex,
    ) -> Vec<(EdgeIndex, T, T)> {
        let eps = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);
        let mut candidates = Vec::new();

        for edge_ref in network.graph.edge_references() {
            let edge_idx = edge_ref.id();
            if edge_idx == previous_edge {
                continue;
            }

            let src = edge_ref.source();
            let dst = edge_ref.target();
            let q = *network.flow_rates.get(&edge_idx).unwrap_or(&T::zero());

            if q.abs() <= eps {
                continue;
            }

            let candidate = if src == node && q > T::zero() {
                Some((edge_idx, T::zero(), q.abs()))
            } else if dst == node && q < T::zero() {
                Some((edge_idx, T::one(), q.abs()))
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
