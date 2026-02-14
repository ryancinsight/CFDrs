use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::solver::{
    EdgeFlowEvent, InletCompositionEvent, MixtureComposition, SimulationTimeConfig,
    PressureBoundaryEvent,
    TransientCompositionSimulator,
};
use std::collections::HashMap;

fn mixture(pairs: &[(i32, f64)]) -> MixtureComposition<f64> {
    let mut map = HashMap::new();
    for (id, frac) in pairs {
        map.insert(*id, *frac);
    }
    MixtureComposition::new(map)
}

#[test]
fn switches_inlet_mixture_over_time() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, 1.0);

    let events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: inlet.index(),
            mixture: mixture(&[(1, 1.0)]),
        },
        InletCompositionEvent {
            time: 2.0,
            node_index: inlet.index(),
            mixture: mixture(&[(2, 1.0)]),
        },
    ];

    let states =
        TransientCompositionSimulator::simulate(&network, events, vec![0.0, 1.0, 2.0, 3.0]).expect("simulate");

    assert_eq!(states.len(), 4);

    let edge0_t1 = states[1].edge_mixtures.get(&edge.index()).expect("edge t1");
    assert!((edge0_t1.fractions.get(&1).copied().unwrap_or(0.0) - 1.0).abs() < 1e-12);

    let edge0_t3 = states[3].edge_mixtures.get(&edge.index()).expect("edge t3");
    assert!((edge0_t3.fractions.get(&2).copied().unwrap_or(0.0) - 1.0).abs() < 1e-12);
}

#[test]
fn mixes_two_inlet_streams_at_junction() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, 2.0);
    network.set_flow_rate(e2, 1.0);
    network.set_flow_rate(e3, 3.0);

    let events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: 0.0,
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let states = TransientCompositionSimulator::simulate(&network, events, vec![0.0]).expect("simulate");
    let out_mix = states[0].edge_mixtures.get(&e3.index()).expect("out edge mixture");

    let f10 = out_mix.fractions.get(&10).copied().unwrap_or(0.0);
    let f20 = out_mix.fractions.get(&20).copied().unwrap_or(0.0);

    assert!((f10 - (2.0 / 3.0)).abs() < 1e-9);
    assert!((f20 - (1.0 / 3.0)).abs() < 1e-9);
}

#[test]
fn timing_config_generates_expected_result_timepoints() {
    let config = SimulationTimeConfig::<f64>::new(1.0, 0.4, 0.1);
    let points = config.result_timepoints().expect("result points");

    assert_eq!(points.len(), 4);
    assert!((points[0] - 0.0).abs() < 1e-12);
    assert!((points[1] - 0.4).abs() < 1e-12);
    assert!((points[2] - 0.8).abs() < 1e-12);
    assert!((points[3] - 1.0).abs() < 1e-12);
}

#[test]
fn simulate_with_time_config_samples_result_grid_and_switches_events() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, 1.0);

    let events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: inlet.index(),
            mixture: mixture(&[(1, 1.0)]),
        },
        InletCompositionEvent {
            time: 0.5,
            node_index: inlet.index(),
            mixture: mixture(&[(2, 1.0)]),
        },
    ];

    let timing = SimulationTimeConfig::new(1.0, 0.25, 0.1);
    let states =
        TransientCompositionSimulator::simulate_with_time_config(&network, events, timing).expect("simulate");

    assert_eq!(states.len(), 5);
    assert!((states[0].time - 0.0).abs() < 1e-12);
    assert!((states[1].time - 0.25).abs() < 1e-12);
    assert!((states[2].time - 0.5).abs() < 1e-12);
    assert!((states[3].time - 0.75).abs() < 1e-12);
    assert!((states[4].time - 1.0).abs() < 1e-12);

    let before_switch = states[1].edge_mixtures.get(&edge.index()).expect("edge t0.25");
    assert!((before_switch.fractions.get(&1).copied().unwrap_or(0.0) - 1.0).abs() < 1e-12);

    let at_switch = states[2].edge_mixtures.get(&edge.index()).expect("edge t0.5");
    assert!((at_switch.fractions.get(&2).copied().unwrap_or(0.0) - 1.0).abs() < 1e-12);
}

#[test]
fn edge_average_concentrations_query_matches_snapshot_mixture() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, 3.0);
    network.set_flow_rate(e2, 1.0);
    network.set_flow_rate(e3, 4.0);

    let events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: 0.0,
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let states = TransientCompositionSimulator::simulate(&network, events, vec![0.0]).expect("simulate");
    let avg = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("average concentrations");

    assert!((avg.get(&10).copied().unwrap_or(0.0) - 0.75).abs() < 1e-9);
    assert!((avg.get(&20).copied().unwrap_or(0.0) - 0.25).abs() < 1e-9);
}

#[test]
fn edge_average_concentrations_query_returns_none_for_missing_edge() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, 1.0);

    let events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: mixture(&[(1, 1.0)]),
    }];

    let states = TransientCompositionSimulator::simulate(&network, events, vec![0.0]).expect("simulate");
    assert!(states[0].average_fluid_concentrations_in_edge(usize::MAX).is_none());
}

#[test]
fn scheduled_flow_events_update_mixture_distribution_over_time() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, 3.0);
    network.set_flow_rate(e2, 1.0);
    network.set_flow_rate(e3, 4.0);

    let composition_events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: 0.0,
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let flow_events = vec![
        EdgeFlowEvent {
            time: 1.0,
            edge_index: e1.index(),
            flow_rate: 1.0,
        },
        EdgeFlowEvent {
            time: 1.0,
            edge_index: e2.index(),
            flow_rate: 3.0,
        },
        EdgeFlowEvent {
            time: 1.0,
            edge_index: e3.index(),
            flow_rate: 4.0,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_flow_events(
        &network,
        composition_events,
        vec![0.0, 1.0],
        flow_events,
    )
    .expect("simulate");

    let before = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("before");
    assert!((before.get(&10).copied().unwrap_or(0.0) - 0.75).abs() < 1e-9);
    assert!((before.get(&20).copied().unwrap_or(0.0) - 0.25).abs() < 1e-9);

    let after = states[1]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("after");
    assert!((after.get(&10).copied().unwrap_or(0.0) - 0.25).abs() < 1e-9);
    assert!((after.get(&20).copied().unwrap_or(0.0) - 0.75).abs() < 1e-9);
}

#[test]
fn scheduled_pressure_events_resolve_flows_and_update_mixture_distribution() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let _e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let _e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("fluid");
    let mut network = Network::new(graph, fluid);

    network.set_pressure(out, 0.0);

    let composition_events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: 0.0,
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: 0.0,
            node_index: in1.index(),
            pressure: 2.5,
        },
        PressureBoundaryEvent {
            time: 0.0,
            node_index: in2.index(),
            pressure: 2.0,
        },
        PressureBoundaryEvent {
            time: 1.0,
            node_index: in1.index(),
            pressure: 2.0,
        },
        PressureBoundaryEvent {
            time: 1.0,
            node_index: in2.index(),
            pressure: 2.5,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        vec![0.0, 1.0],
        pressure_events,
    )
    .expect("simulate");

    let before = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("before concentrations");
    assert!((before.get(&10).copied().unwrap_or(0.0) - (2.0 / 3.0)).abs() < 1e-9);
    assert!((before.get(&20).copied().unwrap_or(0.0) - (1.0 / 3.0)).abs() < 1e-9);

    let after = states[1]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("after concentrations");
    assert!((after.get(&10).copied().unwrap_or(0.0) - (1.0 / 3.0)).abs() < 1e-9);
    assert!((after.get(&20).copied().unwrap_or(0.0) - (2.0 / 3.0)).abs() < 1e-9);
}
