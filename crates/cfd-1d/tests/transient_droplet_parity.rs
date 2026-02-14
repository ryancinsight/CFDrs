use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::solver::{
    CompositionState, DropletInjection, DropletState, EdgeFlowEvent, InletCompositionEvent,
    MixtureComposition, DropletSplitPolicy, PressureBoundaryEvent, SplitMode,
    TransientCompositionSimulator, TransientDropletSimulator,
};
use std::collections::HashMap;

fn pure(fluid_id: i32) -> MixtureComposition<f64> {
    let mut map = HashMap::new();
    map.insert(fluid_id, 1.0);
    MixtureComposition::new(map)
}

fn composition_over_times(
    network: &Network<f64>,
    inlet_node: usize,
    times: Vec<f64>,
) -> Vec<CompositionState<f64>> {
    let events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet_node,
        mixture: pure(1),
    }];
    TransientCompositionSimulator::simulate(network, events, times).expect("composition")
}

#[test]
fn droplet_transitions_to_sink_with_channel_occupancy_tracking() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let out = builder.add_outlet("out".to_string());
    let e0 = builder.connect_with_pipe(inlet, out, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 1.0);

    let composition = composition_over_times(&network, inlet.index(), vec![0.0, 0.4, 1.2]);

    let injections = vec![DropletInjection {
        droplet_id: 7,
        fluid_id: 1,
        volume: 1e-12,
        injection_time: 0.2,
        channel_index: e0.index(),
        relative_position: 0.1,
    }];

    let states =
        TransientDropletSimulator::simulate_on_composition(&network, injections, composition).expect("droplets");

    assert_eq!(states[0].droplets[&7].state, DropletState::Injection);
    assert_eq!(states[1].droplets[&7].state, DropletState::Network);
    assert!(!states[1].droplets[&7].occupied_channels.is_empty());
    assert!(!states[1].droplets[&7].boundaries.is_empty());
    assert!(!states[1].droplets[&7].occupancy_spans.is_empty());
    assert_eq!(states[2].droplets[&7].state, DropletState::Sink);
    assert!(states[2].droplets[&7].occupied_channels.is_empty());
}

#[test]
fn droplet_transitions_to_trapped_when_no_outgoing_path() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let dead = builder.add_junction("dead".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let e0 = builder.connect_with_pipe(inlet, dead, "e0".to_string());
    let e1 = builder.connect_with_pipe(outlet, dead, "e1".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 1.0);
    network.set_flow_rate(e1, 0.5);

    let composition = composition_over_times(&network, inlet.index(), vec![0.0, 0.6, 1.4]);

    let injections = vec![DropletInjection {
        droplet_id: 8,
        fluid_id: 1,
        volume: 1e-12,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.5,
    }];

    let states =
        TransientDropletSimulator::simulate_on_composition(&network, injections, composition).expect("droplets");

    assert_eq!(states[0].droplets[&8].state, DropletState::Network);
    assert_eq!(states[2].droplets[&8].state, DropletState::Trapped);
    assert!(states[2].droplets[&8].position.is_none());
}

#[test]
fn droplet_splits_flow_weighted_with_volume_conservation() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let j = builder.add_junction("j".to_string());
    let out1 = builder.add_outlet("out1".to_string());
    let out2 = builder.add_outlet("out2".to_string());
    let e0 = builder.connect_with_pipe(inlet, j, "e0".to_string());
    let e1 = builder.connect_with_pipe(j, out1, "e1".to_string());
    let e2 = builder.connect_with_pipe(j, out2, "e2".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 2.0);
    network.set_flow_rate(e1, 1.5);
    network.set_flow_rate(e2, 0.5);

    let composition = composition_over_times(&network, inlet.index(), vec![0.0, 0.2, 0.8]);
    let injections = vec![DropletInjection {
        droplet_id: 11,
        fluid_id: 1,
        volume: 0.3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.2,
    }];

    let states =
        TransientDropletSimulator::simulate_on_composition(&network, injections, composition).expect("droplets");

    let s2 = &states[2].droplets[&11];
    assert_eq!(s2.state, DropletState::Network);
    assert!(s2.occupied_channels.contains(&e1.index()));
    assert!(s2.occupied_channels.contains(&e2.index()));
    assert!((s2.total_volume - 0.3).abs() < 1e-9);
}

#[test]
fn split_branches_merge_back_at_reconvergence() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let j0 = builder.add_junction("j0".to_string());
    let j1 = builder.add_junction("j1".to_string());
    let j2 = builder.add_junction("j2".to_string());
    let jm = builder.add_junction("jm".to_string());
    let out = builder.add_outlet("out".to_string());

    let e0 = builder.connect_with_pipe(inlet, j0, "e0".to_string());
    let e1 = builder.connect_with_pipe(j0, j1, "e1".to_string());
    let e2 = builder.connect_with_pipe(j0, j2, "e2".to_string());
    let e3 = builder.connect_with_pipe(j1, jm, "e3".to_string());
    let e4 = builder.connect_with_pipe(j2, jm, "e4".to_string());
    let e5 = builder.connect_with_pipe(jm, out, "e5".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 2.0);
    network.set_flow_rate(e1, 1.0);
    network.set_flow_rate(e2, 1.0);
    network.set_flow_rate(e3, 1.0);
    network.set_flow_rate(e4, 1.0);
    network.set_flow_rate(e5, 2.0);

    let composition = composition_over_times(&network, inlet.index(), vec![0.0, 0.6, 1.6, 2.6]);
    let injections = vec![DropletInjection {
        droplet_id: 12,
        fluid_id: 1,
        volume: 0.2,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.1,
    }];

    let states =
        TransientDropletSimulator::simulate_on_composition(&network, injections, composition).expect("droplets");

    let split_state = &states[1].droplets[&12];
    assert_eq!(split_state.state, DropletState::Network);
    assert!(split_state.occupied_channels.len() >= 2);

    let merged_state = &states[3].droplets[&12];
    assert_eq!(merged_state.state, DropletState::Network);
    assert_eq!(merged_state.occupied_channels.len(), 1);
    assert_eq!(merged_state.occupied_channels[0], e5.index());
    assert!((merged_state.total_volume - 0.2).abs() < 1e-9);
}

#[test]
fn auto_policy_uses_no_split_for_dominant_branch_scenarios() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let j = builder.add_junction("j".to_string());
    let out1 = builder.add_outlet("out1".to_string());
    let out2 = builder.add_outlet("out2".to_string());
    let e0 = builder.connect_with_pipe(inlet, j, "e0".to_string());
    let e1 = builder.connect_with_pipe(j, out1, "e1".to_string());
    let e2 = builder.connect_with_pipe(j, out2, "e2".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 2.0);
    network.set_flow_rate(e1, 1.9);
    network.set_flow_rate(e2, 0.1);

    let composition = composition_over_times(&network, inlet.index(), vec![0.0, 0.2, 0.8]);
    let injections = vec![DropletInjection {
        droplet_id: 21,
        fluid_id: 1,
        volume: 0.3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.2,
    }];

    let policy = DropletSplitPolicy {
        mode: SplitMode::AutoFlowWeighted,
        min_secondary_flow_fraction: 0.2,
        min_child_volume: 1e-15,
        max_split_branches: 2,
    };

    let states = TransientDropletSimulator::simulate_on_composition_with_policy(
        &network,
        injections,
        composition,
        policy,
    )
    .expect("droplets");

    let s2 = &states[2].droplets[&21];
    assert_eq!(s2.state, DropletState::Network);
    assert_eq!(s2.occupied_channels.len(), 1);
    assert_eq!(s2.occupied_channels[0], e1.index());
    assert!((s2.total_volume - 0.3).abs() < 1e-9);
}

#[test]
fn pressure_event_droplet_api_matches_manual_composition_pipeline() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let out = builder.add_outlet("out".to_string());
    let e0 = builder.connect_with_pipe(inlet, out, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, 0.0);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: 0.0,
            node_index: inlet.index(),
            pressure: 1.0,
        },
        PressureBoundaryEvent {
            time: 0.5,
            node_index: inlet.index(),
            pressure: 2.0,
        },
    ];

    let timepoints = vec![0.0, 0.2, 0.6];
    let injections = vec![DropletInjection {
        droplet_id: 40,
        fluid_id: 1,
        volume: 1e-3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.1,
    }];

    let states_api = TransientDropletSimulator::simulate_with_pressure_events(
        &network,
        injections.clone(),
        composition_events.clone(),
        timepoints.clone(),
        pressure_events.clone(),
    )
    .expect("pressure-event droplet api");

    let composition_states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        timepoints,
        pressure_events,
    )
    .expect("manual composition");
    let states_manual = TransientDropletSimulator::simulate_on_composition(
        &network,
        injections,
        composition_states,
    )
    .expect("manual droplet");

    assert_eq!(states_api.len(), states_manual.len());
    for (api, manual) in states_api.iter().zip(states_manual.iter()) {
        assert!((api.time - manual.time).abs() < 1e-12);
        let d_api = &api.droplets[&40];
        let d_manual = &manual.droplets[&40];
        assert_eq!(d_api.state, d_manual.state);
        assert_eq!(d_api.occupied_channels, d_manual.occupied_channels);
        assert!((d_api.total_volume - d_manual.total_volume).abs() < 1e-12);
    }
}

#[test]
fn higher_pressure_event_drives_faster_droplet_progress() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let out = builder.add_outlet("out".to_string());
    let e0 = builder.connect_with_pipe(inlet, out, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, 0.0);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];
    let timepoints = vec![0.0, 0.2];

    let injections = vec![DropletInjection {
        droplet_id: 41,
        fluid_id: 1,
        volume: 1e-3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.1,
    }];

    let low_pressure_states = TransientDropletSimulator::simulate_with_pressure_events(
        &network,
        injections.clone(),
        composition_events.clone(),
        timepoints.clone(),
        vec![PressureBoundaryEvent {
            time: 0.0,
            node_index: inlet.index(),
            pressure: 1.0,
        }],
    )
    .expect("low pressure run");

    let high_pressure_states = TransientDropletSimulator::simulate_with_pressure_events(
        &network,
        injections,
        composition_events,
        timepoints,
        vec![PressureBoundaryEvent {
            time: 0.0,
            node_index: inlet.index(),
            pressure: 2.0,
        }],
    )
    .expect("high pressure run");

    let low_pos = low_pressure_states[1].droplets[&41]
        .position
        .as_ref()
        .expect("low pressure position")
        .relative_position;
    let high_pos = high_pressure_states[1].droplets[&41]
        .position
        .as_ref()
        .expect("high pressure position")
        .relative_position;

    assert!(high_pos > low_pos);
}

#[test]
fn flow_event_droplet_api_matches_manual_composition_pipeline() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let out = builder.add_outlet("out".to_string());
    let e0 = builder.connect_with_pipe(inlet, out, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 1.0);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];
    let flow_events = vec![
        EdgeFlowEvent {
            time: 0.0,
            edge_index: e0.index(),
            flow_rate: 1.0,
        },
        EdgeFlowEvent {
            time: 0.5,
            edge_index: e0.index(),
            flow_rate: 2.0,
        },
    ];
    let timepoints = vec![0.0, 0.2, 0.6];

    let injections = vec![DropletInjection {
        droplet_id: 50,
        fluid_id: 1,
        volume: 1e-3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.1,
    }];

    let states_api = TransientDropletSimulator::simulate_with_flow_events(
        &network,
        injections.clone(),
        composition_events.clone(),
        timepoints.clone(),
        flow_events.clone(),
    )
    .expect("flow-event droplet api");

    let composition_states = TransientCompositionSimulator::simulate_with_flow_events(
        &network,
        composition_events,
        timepoints,
        flow_events,
    )
    .expect("manual composition");
    let states_manual = TransientDropletSimulator::simulate_on_composition(
        &network,
        injections,
        composition_states,
    )
    .expect("manual droplet");

    assert_eq!(states_api.len(), states_manual.len());
    for (api, manual) in states_api.iter().zip(states_manual.iter()) {
        assert!((api.time - manual.time).abs() < 1e-12);
        let d_api = &api.droplets[&50];
        let d_manual = &manual.droplets[&50];
        assert_eq!(d_api.state, d_manual.state);
        assert_eq!(d_api.occupied_channels, d_manual.occupied_channels);
        assert!((d_api.total_volume - d_manual.total_volume).abs() < 1e-12);
    }
}

#[test]
fn end_to_end_flow_event_policy_controls_branching_behavior() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let j = builder.add_junction("j".to_string());
    let out1 = builder.add_outlet("out1".to_string());
    let out2 = builder.add_outlet("out2".to_string());
    let e0 = builder.connect_with_pipe(inlet, j, "e0".to_string());
    let e1 = builder.connect_with_pipe(j, out1, "e1".to_string());
    let e2 = builder.connect_with_pipe(j, out2, "e2".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e0, 2.0);
    network.set_flow_rate(e1, 1.9);
    network.set_flow_rate(e2, 0.1);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];
    let flow_events = vec![
        EdgeFlowEvent {
            time: 0.3,
            edge_index: e1.index(),
            flow_rate: 1.0,
        },
        EdgeFlowEvent {
            time: 0.3,
            edge_index: e2.index(),
            flow_rate: 1.0,
        },
        EdgeFlowEvent {
            time: 0.3,
            edge_index: e0.index(),
            flow_rate: 2.0,
        },
    ];

    let injections = vec![DropletInjection {
        droplet_id: 51,
        fluid_id: 1,
        volume: 0.3,
        injection_time: 0.0,
        channel_index: e0.index(),
        relative_position: 0.2,
    }];

    let policy = DropletSplitPolicy {
        mode: SplitMode::AutoFlowWeighted,
        min_secondary_flow_fraction: 0.2,
        min_child_volume: 1e-15,
        max_split_branches: 2,
    };

    let states = TransientDropletSimulator::simulate_with_flow_events_and_policy(
        &network,
        injections,
        composition_events,
        vec![0.0, 0.2, 0.6],
        flow_events,
        policy,
    )
    .expect("flow event droplet run");

    let pre_switch = &states[1].droplets[&51];
    assert_eq!(pre_switch.state, DropletState::Network);
    assert_eq!(pre_switch.occupied_channels.len(), 1);
    assert_eq!(pre_switch.occupied_channels[0], e0.index());

    let post_switch = &states[2].droplets[&51];
    assert_eq!(post_switch.state, DropletState::Network);
    assert!(post_switch.occupied_channels.contains(&e1.index()));
    assert!(post_switch.occupied_channels.contains(&e2.index()));
    assert!((post_switch.total_volume - 0.3).abs() < 1e-9);
}
