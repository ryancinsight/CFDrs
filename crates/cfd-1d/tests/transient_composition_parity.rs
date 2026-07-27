use aequitas::systems::si::quantities::{
    Area, Dimensionless, HydraulicResistance, Length, Pressure, Time, VolumetricFlowRate,
};
use cfd_1d::domain::network::{EdgeProperties, Network, NetworkBuilder, ResistanceUpdatePolicy};
use cfd_1d::pries_phase_separation;
use cfd_1d::solver::core::{
    BloodEdgeTransportConfig, EdgeFlowEvent, InletCompositionEvent, InletHematocritEvent,
    MixtureComposition, PressureBoundaryEvent, SimulationTimeConfig, TransientCompositionSimulator,
};
use cfd_1d::{ChannelGeometry, CrossSection, SurfaceProperties, Wettability};
use std::collections::HashMap;

fn mixture(pairs: &[(i32, f64)]) -> MixtureComposition<f64> {
    let mut map = HashMap::new();
    for (id, frac) in pairs {
        map.insert(*id, fraction(*frac));
    }
    MixtureComposition::new(map)
}

fn fraction(value: f64) -> Dimensionless<f64> {
    Dimensionless::from_base(value)
}

fn concentration(map: &HashMap<i32, Dimensionless<f64>>, fluid_id: i32) -> f64 {
    map.get(&fluid_id)
        .copied()
        .map(Dimensionless::into_base)
        .unwrap_or(0.0)
}

fn time(value: f64) -> Time<f64> {
    Time::from_base(value)
}

fn timepoints(values: &[f64]) -> Vec<Time<f64>> {
    values.iter().copied().map(time).collect()
}

fn flow(value: f64) -> VolumetricFlowRate<f64> {
    VolumetricFlowRate::from_base(value)
}

fn pressure(value: f64) -> Pressure<f64> {
    Pressure::from_base(value)
}

fn cstr_relaxation(previous: f64, inlet: f64, dt: f64, tau: f64) -> f64 {
    inlet + (previous - inlet) * (-(dt / tau)).exp()
}

fn segmented_transport_config() -> BloodEdgeTransportConfig<f64> {
    BloodEdgeTransportConfig::new(4, fraction(1.0))
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
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));

    let events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: inlet.index(),
            mixture: mixture(&[(1, 1.0)]),
        },
        InletCompositionEvent {
            time: time(2.0),
            node_index: inlet.index(),
            mixture: mixture(&[(2, 1.0)]),
        },
    ];

    let states = TransientCompositionSimulator::simulate(
        &network,
        events,
        timepoints(&[0.0, 1.0, 2.0, 3.0]),
    )
    .expect("simulate");

    assert_eq!(states.len(), 4);

    let edge0_t1 = states[1].edge_mixtures.get(&edge.index()).expect("edge t1");
    assert!((concentration(&edge0_t1.fractions, 1) - 1.0).abs() < 1e-12);

    let edge0_t3 = states[3].edge_mixtures.get(&edge.index()).expect("edge t3");
    assert!((concentration(&edge0_t3.fractions, 2) - 1.0).abs() < 1e-12);
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
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(2.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(3.0));

    let events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: time(0.0),
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let states = TransientCompositionSimulator::simulate(&network, events, timepoints(&[0.0]))
        .expect("simulate");
    let out_mix = states[0]
        .edge_mixtures
        .get(&e3.index())
        .expect("out edge mixture");

    let f10 = concentration(&out_mix.fractions, 10);
    let f20 = concentration(&out_mix.fractions, 20);

    assert!((f10 - (2.0 / 3.0)).abs() < 1e-9);
    assert!((f20 - (1.0 / 3.0)).abs() < 1e-9);
}

#[test]
fn timing_config_generates_expected_result_timepoints() {
    let config = SimulationTimeConfig::<f64>::new(time(1.0), time(0.4), time(0.1));
    let points = config.result_timepoints().expect("result points");

    assert_eq!(points.len(), 4);
    assert!((points[0].into_base() - 0.0).abs() < 1e-12);
    assert!((points[1].into_base() - 0.4).abs() < 1e-12);
    assert!((points[2].into_base() - 0.8).abs() < 1e-12);
    assert!((points[3].into_base() - 1.0).abs() < 1e-12);
}

#[test]
fn typed_composition_control_metrics_preserve_values() {
    let timing = SimulationTimeConfig::new(time(2.0), time(0.5), time(0.1));
    let hematocrit_event = InletHematocritEvent {
        time: time(0.25),
        node_index: 0,
        hematocrit: fraction(0.45),
    };
    let flow_event = EdgeFlowEvent {
        time: time(0.5),
        edge_index: 1,
        flow_rate: flow(2.5),
    };
    let pressure_event = PressureBoundaryEvent {
        time: time(0.75),
        node_index: 0,
        pressure: pressure(12.0),
    };
    let transport = BloodEdgeTransportConfig::new(4, fraction(0.8));

    assert_eq!(timing.duration.into_base(), 2.0);
    assert_eq!(hematocrit_event.hematocrit.into_base(), 0.45);
    assert_eq!(flow_event.flow_rate.into_base(), 2.5);
    assert_eq!(pressure_event.pressure.into_base(), 12.0);
    assert_eq!(transport.max_courant_number.into_base(), 0.8);
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
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));

    let events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: inlet.index(),
            mixture: mixture(&[(1, 1.0)]),
        },
        InletCompositionEvent {
            time: time(0.5),
            node_index: inlet.index(),
            mixture: mixture(&[(2, 1.0)]),
        },
    ];

    let timing = SimulationTimeConfig::new(time(1.0), time(0.25), time(0.1));
    let states = TransientCompositionSimulator::simulate_with_time_config(&network, events, timing)
        .expect("simulate");

    assert_eq!(states.len(), 5);
    assert!((states[0].time.into_base() - 0.0).abs() < 1e-12);
    assert!((states[1].time.into_base() - 0.25).abs() < 1e-12);
    assert!((states[2].time.into_base() - 0.5).abs() < 1e-12);
    assert!((states[3].time.into_base() - 0.75).abs() < 1e-12);
    assert!((states[4].time.into_base() - 1.0).abs() < 1e-12);

    let before_switch = states[1]
        .edge_mixtures
        .get(&edge.index())
        .expect("edge t0.25");
    assert!((concentration(&before_switch.fractions, 1) - 1.0).abs() < 1e-12);

    let at_switch = states[2]
        .edge_mixtures
        .get(&edge.index())
        .expect("edge t0.5");
    assert!((concentration(&at_switch.fractions, 2) - 1.0).abs() < 1e-12);
}

#[test]
fn duplicate_requested_timepoints_are_sampled_consistently() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));

    let events = vec![InletHematocritEvent {
        time: time(0.0),
        node_index: inlet.index(),
        hematocrit: fraction(0.45),
    }];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit_with_edge_transport(
        &network,
        events,
        timepoints(&[0.0, 0.5, 0.5, 1.0]),
    )
    .expect("simulate duplicate timepoints");

    assert_eq!(states.len(), 4);
    assert!((states[1].time.into_base() - 0.5).abs() < 1.0e-12);
    assert!((states[2].time.into_base() - 0.5).abs() < 1.0e-12);

    let hct_1 = states[1]
        .edge_hematocrit(edge.index())
        .map(Dimensionless::into_base)
        .expect("edge hct");
    let hct_2 = states[2]
        .edge_hematocrit(edge.index())
        .map(Dimensionless::into_base)
        .expect("edge hct");
    assert!((hct_1 - hct_2).abs() < 1.0e-12);
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
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(3.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(4.0));

    let events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: time(0.0),
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let states = TransientCompositionSimulator::simulate(&network, events, timepoints(&[0.0]))
        .expect("simulate");
    let avg = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("average concentrations");

    assert!((concentration(&avg, 10) - 0.75).abs() < 1e-9);
    assert!((concentration(&avg, 20) - 0.25).abs() < 1e-9);
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
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));

    let events = vec![InletCompositionEvent {
        time: time(0.0),
        node_index: inlet.index(),
        mixture: mixture(&[(1, 1.0)]),
    }];

    let states = TransientCompositionSimulator::simulate(&network, events, timepoints(&[0.0]))
        .expect("simulate");
    assert!(states[0]
        .average_fluid_concentrations_in_edge(usize::MAX)
        .is_none());
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
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(3.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(4.0));

    let composition_events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: time(0.0),
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let flow_events = vec![
        EdgeFlowEvent {
            time: time(1.0),
            edge_index: e1.index(),
            flow_rate: flow(1.0),
        },
        EdgeFlowEvent {
            time: time(1.0),
            edge_index: e2.index(),
            flow_rate: flow(3.0),
        },
        EdgeFlowEvent {
            time: time(1.0),
            edge_index: e3.index(),
            flow_rate: flow(4.0),
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_flow_events(
        &network,
        composition_events,
        timepoints(&[0.0, 1.0]),
        flow_events,
    )
    .expect("simulate");

    let before = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("before");
    assert!((concentration(&before, 10) - 0.75).abs() < 1e-9);
    assert!((concentration(&before, 20) - 0.25).abs() < 1e-9);

    let after = states[1]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("after");
    assert!((concentration(&after, 10) - 0.25).abs() < 1e-9);
    assert!((concentration(&after, 20) - 0.75).abs() < 1e-9);
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

    network.set_pressure(out, Pressure::from_base(0.0));

    let composition_events = vec![
        InletCompositionEvent {
            time: time(0.0),
            node_index: in1.index(),
            mixture: mixture(&[(10, 1.0)]),
        },
        InletCompositionEvent {
            time: time(0.0),
            node_index: in2.index(),
            mixture: mixture(&[(20, 1.0)]),
        },
    ];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in1.index(),
            pressure: pressure(2.5),
        },
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in2.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in1.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in2.index(),
            pressure: pressure(2.5),
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        timepoints(&[0.0, 1.0]),
        pressure_events,
    )
    .expect("simulate");

    let before = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("before concentrations");
    assert!((concentration(&before, 10) - (2.0 / 3.0)).abs() < 1e-9);
    assert!((concentration(&before, 20) - (1.0 / 3.0)).abs() < 1e-9);

    let after = states[1]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("after concentrations");
    assert!((concentration(&after, 10) - (1.0 / 3.0)).abs() < 1e-9);
    assert!((concentration(&after, 20) - (2.0 / 3.0)).abs() < 1e-9);
}

#[test]
fn blood_hematocrit_switches_over_time() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));

    let events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: inlet.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(2.0),
            node_index: inlet.index(),
            hematocrit: fraction(0.20),
        },
    ];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit(
        &network,
        events,
        timepoints(&[0.0, 1.0, 2.0, 3.0]),
    )
    .expect("simulate blood hematocrit");

    assert!(
        (states[1]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.45)
            .abs()
            < 1e-12
    );
    assert!(
        (states[3]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.20)
            .abs()
            < 1e-12
    );
}

#[test]
fn blood_hematocrit_mixes_two_inlets_at_junction() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(3.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(4.0));

    let events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.15),
        },
    ];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit(
        &network,
        events,
        timepoints(&[0.0]),
    )
    .expect("simulate blood hematocrit");
    let out_hct = states[0]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .expect("outlet hematocrit");

    assert!((out_hct - 0.375).abs() < 1e-9);
}

#[test]
fn pressure_event_blood_hematocrit_tracks_flow_redistribution() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let _e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let _e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, Pressure::from_base(0.0));

    let hematocrit_events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.15),
        },
    ];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in1.index(),
            pressure: pressure(2.5),
        },
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in2.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in1.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in2.index(),
            pressure: pressure(2.5),
        },
    ];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit_with_pressure_events(
        &network,
        hematocrit_events,
        timepoints(&[0.0, 1.0]),
        pressure_events,
    )
    .expect("simulate blood hematocrit");

    let before = states[0]
        .edge_hematocrit(e3.index())
        .expect("before hematocrit");
    let after = states[1]
        .edge_hematocrit(e3.index())
        .expect("after hematocrit");

    assert!((before.into_base() - 0.35).abs() < 1e-9);
    assert!((after.into_base() - 0.25).abs() < 1e-9);
}

#[test]
fn coupled_pressure_event_blood_hematocrit_feeds_back_on_flow_split() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, Pressure::from_base(0.0));

    let inlet_diameter = 70.0e-6;
    let inlet_length = 20.0e-3;
    let outlet_diameter = 140.0e-6;
    let outlet_length = 2.0e-3;
    let make_props = |id: &str, diameter: f64, length: f64| EdgeProperties {
        id: id.to_string(),
        component_type: cfd_1d::ComponentType::Pipe,
        length: Length::from_base(length),
        area: Area::from_base(std::f64::consts::PI * diameter * diameter / 4.0),
        hydraulic_diameter: Some(Length::from_base(diameter)),
        resistance: HydraulicResistance::from_base(0.0),
        geometry: Some(ChannelGeometry {
            channel_type: cfd_1d::ChannelType::Straight,
            length: Length::from_base(length),
            cross_section: CrossSection::Circular {
                diameter: Length::from_base(diameter),
            },
            surface: SurfaceProperties {
                roughness: Length::from_base(0.0),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }),
        resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
        properties: HashMap::new(),
    };
    network.add_edge_properties(e1, make_props("e1", inlet_diameter, inlet_length));
    network.add_edge_properties(e2, make_props("e2", inlet_diameter, inlet_length));
    network.add_edge_properties(e3, make_props("e3", outlet_diameter, outlet_length));

    let hematocrit_events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.10),
        },
    ];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in1.index(),
            pressure: pressure(0.02),
        },
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in2.index(),
            pressure: pressure(0.02),
        },
    ];

    let uncoupled = TransientCompositionSimulator::simulate_blood_hematocrit_with_pressure_events(
        &network,
        hematocrit_events.clone(),
        timepoints(&[0.0]),
        pressure_events.clone(),
    )
    .expect("uncoupled blood simulation");
    let coupled =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_coupled_pressure_events(
            &network,
            hematocrit_events,
            timepoints(&[0.0]),
            pressure_events,
        )
        .expect("coupled blood simulation");

    let coupled_hct = coupled[0]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .expect("coupled outlet hematocrit");
    let uncoupled_q1 = uncoupled[0]
        .edge_flow_rates
        .get(&e1.index())
        .copied()
        .map(|q| q.into_base())
        .expect("uncoupled branch 1 flow");
    let uncoupled_q2 = uncoupled[0]
        .edge_flow_rates
        .get(&e2.index())
        .copied()
        .map(|q| q.into_base())
        .expect("uncoupled branch 2 flow");
    let coupled_q1 = coupled[0]
        .edge_flow_rates
        .get(&e1.index())
        .copied()
        .map(|q| q.into_base())
        .expect("coupled branch 1 flow");
    let coupled_q2 = coupled[0]
        .edge_flow_rates
        .get(&e2.index())
        .copied()
        .map(|q| q.into_base())
        .expect("coupled branch 2 flow");
    let coupled_h1 = coupled[0]
        .edge_hematocrit(e1.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);
    let coupled_h2 = coupled[0]
        .edge_hematocrit(e2.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);

    assert!((uncoupled_q1 - uncoupled_q2).abs() < 1.0e-12);
    assert!(coupled_q1.is_finite() && coupled_q1 > 0.0);
    assert!(coupled_q2.is_finite() && coupled_q2 > 0.0);
    assert!((coupled_h1 - 0.45).abs() < 1.0e-12);
    assert!((coupled_h2 - 0.10).abs() < 1.0e-12);
    assert!(coupled_hct.is_finite());
    assert!((0.10..=0.45).contains(&coupled_hct));
}

#[test]
fn blood_edge_transport_relaxes_single_channel_hematocrit_front() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));
    network.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e0".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: Length::from_base(1.0),
            area: Area::from_base(1.0),
            hydraulic_diameter: Some(Length::from_base(1.0)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: Length::from_base(1.0),
                cross_section: CrossSection::Circular {
                    diameter: Length::from_base(1.0),
                },
                surface: SurfaceProperties {
                    roughness: Length::from_base(0.0),
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        },
    );

    let events = vec![InletHematocritEvent {
        time: time(0.0),
        node_index: inlet.index(),
        hematocrit: fraction(0.45),
    }];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit_with_edge_transport(
        &network,
        events,
        timepoints(&[0.0, 0.5, 1.0, 2.0]),
    )
    .expect("simulate transported blood hematocrit");

    let expected_t05 = cstr_relaxation(0.0, 0.45, 0.5, 1.0);
    let expected_t10 = cstr_relaxation(expected_t05, 0.45, 0.5, 1.0);
    let expected_t20 = cstr_relaxation(expected_t10, 0.45, 1.0, 1.0);

    assert!(
        states[0]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(-1.0)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[1]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t05)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[2]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t10)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[3]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t20)
            .abs()
            < 1.0e-12
    );
}

#[test]
fn blood_edge_transport_relaxes_mixed_outlet_hematocrit() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(3.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(4.0));
    network.add_edge_properties(
        e3,
        EdgeProperties {
            id: "e3".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: Length::from_base(2.0),
            area: Area::from_base(2.0),
            hydraulic_diameter: Some(Length::from_base(1.0)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: Length::from_base(2.0),
                cross_section: CrossSection::Circular {
                    diameter: Length::from_base(1.0),
                },
                surface: SurfaceProperties {
                    roughness: Length::from_base(0.0),
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        },
    );

    let events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.15),
        },
    ];

    let states = TransientCompositionSimulator::simulate_blood_hematocrit_with_edge_transport(
        &network,
        events,
        timepoints(&[0.0, 0.5, 1.0]),
    )
    .expect("simulate transported mixed blood hematocrit");

    let expected_t05 = cstr_relaxation(0.0, 0.375, 0.5, 1.0);
    let expected_t10 = cstr_relaxation(expected_t05, 0.375, 0.5, 1.0);

    assert!(
        states[0]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(-1.0)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[1]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t05)
            .abs()
            < 1.0e-9
    );
    assert!(
        (states[2]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t10)
            .abs()
            < 1.0e-9
    );
}

#[test]
fn blood_pressure_event_edge_transport_uses_edge_inventory() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let _e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let _e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, Pressure::from_base(0.0));
    network.add_edge_properties(
        e3,
        EdgeProperties {
            id: "e3".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: Length::from_base(1.0),
            area: Area::from_base(1.0),
            hydraulic_diameter: Some(Length::from_base(1.0)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: Length::from_base(1.0),
                cross_section: CrossSection::Circular {
                    diameter: Length::from_base(1.0),
                },
                surface: SurfaceProperties {
                    roughness: Length::from_base(0.0),
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        },
    );

    let hematocrit_events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.15),
        },
    ];
    let pressure_events = vec![
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in1.index(),
            pressure: pressure(2.5),
        },
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in2.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in1.index(),
            pressure: pressure(2.0),
        },
        PressureBoundaryEvent {
            time: time(1.0),
            node_index: in2.index(),
            pressure: pressure(2.5),
        },
    ];

    let instant_states =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_pressure_events(
            &network,
            hematocrit_events.clone(),
            timepoints(&[0.0, 1.0, 2.0]),
            pressure_events.clone(),
        )
        .expect("simulate instantaneous pressure-event blood transport");
    let states =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_pressure_events_and_edge_transport(
            &network,
            hematocrit_events,
            timepoints(&[0.0, 1.0, 2.0]),
            pressure_events,
        )
        .expect("simulate pressure-event edge transport");

    let instant_t0 = instant_states[0]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .unwrap_or(0.0);
    let instant_t1 = instant_states[1]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .unwrap_or(0.0);
    let tau_t0 = 1.0
        / instant_states[0]
            .edge_flow_rates
            .get(&e3.index())
            .copied()
            .map(|q| q.into_base())
            .unwrap_or(1.0)
            .abs();
    let tau_t1 = 1.0
        / instant_states[1]
            .edge_flow_rates
            .get(&e3.index())
            .copied()
            .map(|q| q.into_base())
            .unwrap_or(1.0)
            .abs();
    let expected_t10 = cstr_relaxation(0.0, instant_t0, 1.0, tau_t0);
    let expected_t20 = cstr_relaxation(expected_t10, instant_t1, 1.0, tau_t1);

    assert!(
        states[0]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(-1.0)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[1]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t10)
            .abs()
            < 1.0e-9
    );
    assert!(
        (states[2]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - expected_t20)
            .abs()
            < 1.0e-9
    );
}

#[test]
fn blood_segmented_edge_transport_advects_front_one_cell_per_step() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(edge, VolumetricFlowRate::from_base(1.0));
    network.add_edge_properties(
        edge,
        EdgeProperties {
            id: "e0".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: Length::from_base(1.0),
            area: Area::from_base(1.0),
            hydraulic_diameter: Some(Length::from_base(1.0)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: Length::from_base(1.0),
                cross_section: CrossSection::Circular {
                    diameter: Length::from_base(1.0),
                },
                surface: SurfaceProperties {
                    roughness: Length::from_base(0.0),
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        },
    );

    let states =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_segmented_edge_transport(
            &network,
            vec![InletHematocritEvent {
                time: time(0.0),
                node_index: inlet.index(),
                hematocrit: fraction(0.45),
            }],
            timepoints(&[0.0, 0.25, 0.50, 0.75, 1.0]),
            segmented_transport_config(),
        )
        .expect("simulate segmented blood transport");

    assert!(
        states[0]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(-1.0)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[1]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.1125)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[2]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.2250)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[3]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.3375)
            .abs()
            < 1.0e-12
    );
    assert!(
        (states[4]
            .edge_hematocrit(edge.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.4500)
            .abs()
            < 1.0e-12
    );
}

#[test]
fn blood_segmented_edge_transport_delays_mixed_outlet_more_than_single_volume_model() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(e1, VolumetricFlowRate::from_base(3.0));
    network.set_flow_rate(e2, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(e3, VolumetricFlowRate::from_base(4.0));
    network.add_edge_properties(
        e3,
        EdgeProperties {
            id: "e3".to_string(),
            component_type: cfd_1d::ComponentType::Pipe,
            length: Length::from_base(2.0),
            area: Area::from_base(2.0),
            hydraulic_diameter: Some(Length::from_base(1.0)),
            resistance: HydraulicResistance::from_base(1.0),
            geometry: Some(ChannelGeometry {
                channel_type: cfd_1d::ChannelType::Straight,
                length: Length::from_base(2.0),
                cross_section: CrossSection::Circular {
                    diameter: Length::from_base(1.0),
                },
                surface: SurfaceProperties {
                    roughness: Length::from_base(0.0),
                    contact_angle: None,
                    surface_energy: None,
                    wettability: Wettability::Hydrophilic,
                },
                variations: Vec::new(),
            }),
            resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
            properties: HashMap::new(),
        },
    );

    let events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.15),
        },
    ];

    let lumped = TransientCompositionSimulator::simulate_blood_hematocrit_with_edge_transport(
        &network,
        events.clone(),
        timepoints(&[0.0, 0.25, 0.50, 0.75, 1.0]),
    )
    .expect("simulate lumped edge transport");
    let segmented =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_segmented_edge_transport(
            &network,
            events,
            timepoints(&[0.0, 0.25, 0.50, 0.75, 1.0]),
            segmented_transport_config(),
        )
        .expect("simulate segmented edge transport");

    assert!(
        (segmented[1]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.09375)
            .abs()
            < 1.0e-12
    );
    assert!(
        (segmented[4]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            - 0.37500)
            .abs()
            < 1.0e-12
    );
    assert!(
        segmented[1]
            .edge_hematocrit(e3.index())
            .map(Dimensionless::into_base)
            .unwrap_or(0.0)
            > lumped[1]
                .edge_hematocrit(e3.index())
                .map(Dimensionless::into_base)
                .unwrap_or(0.0)
    );
}

#[test]
fn blood_segmented_edge_transport_applies_pries_split_at_bifurcation() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let junction = builder.add_junction("junction".to_string());
    let small_outlet = builder.add_outlet("small_outlet".to_string());
    let large_outlet = builder.add_outlet("large_outlet".to_string());

    let parent = builder.connect_with_pipe(inlet, junction, "parent".to_string());
    let small = builder.connect_with_pipe(junction, small_outlet, "small".to_string());
    let large = builder.connect_with_pipe(junction, large_outlet, "large".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_flow_rate(parent, VolumetricFlowRate::from_base(1.0));
    network.set_flow_rate(small, VolumetricFlowRate::from_base(0.5));
    network.set_flow_rate(large, VolumetricFlowRate::from_base(0.5));

    let make_props = |id: &str, hydraulic_diameter: f64| EdgeProperties {
        id: id.to_string(),
        component_type: cfd_1d::ComponentType::Pipe,
        length: Length::from_base(1.0),
        area: Area::from_base(1.0),
        hydraulic_diameter: Some(Length::from_base(hydraulic_diameter)),
        resistance: HydraulicResistance::from_base(1.0),
        geometry: Some(ChannelGeometry {
            channel_type: cfd_1d::ChannelType::Straight,
            length: Length::from_base(1.0),
            cross_section: CrossSection::Circular {
                diameter: Length::from_base(hydraulic_diameter),
            },
            surface: SurfaceProperties {
                roughness: Length::from_base(0.0),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }),
        resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
        properties: HashMap::new(),
    };

    network.add_edge_properties(parent, make_props("parent", 60.0e-6));
    network.add_edge_properties(small, make_props("small", 40.0e-6));
    network.add_edge_properties(large, make_props("large", 80.0e-6));
    network
        .properties
        .get_mut(&parent)
        .expect("parent props")
        .properties
        .insert(cfd_1d::EDGE_PROPERTY_LOCAL_HEMATOCRIT.to_string(), 0.45);

    let states =
        TransientCompositionSimulator::simulate_blood_hematocrit_with_segmented_edge_transport(
            &network,
            vec![InletHematocritEvent {
                time: time(0.0),
                node_index: inlet.index(),
                hematocrit: fraction(0.45),
            }],
            timepoints(&[0.0, 0.25]),
            BloodEdgeTransportConfig::new(1, fraction(1.0)),
        )
        .expect("simulate bifurcation segmented transport");

    let expected_small = pries_phase_separation(
        0.45,
        0.5,
        Length::from_base(40.0e-6),
        Length::from_base(80.0e-6),
        Length::from_base(60.0e-6),
    )
    .expect("pries split");
    let expected_large_inlet = ((1.0 - expected_small.cell_fraction) * 0.45 / 0.5).clamp(0.0, 1.0);
    let step_factor = 0.25 / (1.0 / 0.5);

    let small_hct = states[1]
        .edge_hematocrit(small.index())
        .map(Dimensionless::into_base)
        .expect("small hct");
    let large_hct = states[1]
        .edge_hematocrit(large.index())
        .map(Dimensionless::into_base)
        .expect("large hct");

    assert!((small_hct - expected_small.daughter_hematocrit * step_factor).abs() < 1.0e-12);
    assert!((large_hct - expected_large_inlet * step_factor).abs() < 1.0e-12);
    assert!((small_hct - large_hct).abs() > 1.0e-12);
    assert!(
        ((small_hct / step_factor) * 0.5 + (large_hct / step_factor) * 0.5 - 0.45).abs() < 1.0e-12
    );
}

#[test]
fn coupled_pressure_event_segmented_blood_transport_feeds_back_on_resistance_updates() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::blood::CarreauYasudaBlood::<f64>::normal_blood();
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, Pressure::from_base(0.0));

    let inlet_diameter = 1.0;
    let inlet_length = 1.0;
    let outlet_diameter = 1.0;
    let outlet_length = 1.0;
    let make_props = |id: &str, diameter: f64, length: f64| EdgeProperties {
        id: id.to_string(),
        component_type: cfd_1d::ComponentType::Pipe,
        length: Length::from_base(length),
        area: Area::from_base(std::f64::consts::PI * diameter * diameter / 4.0),
        hydraulic_diameter: Some(Length::from_base(diameter)),
        resistance: HydraulicResistance::from_base(0.0),
        geometry: Some(ChannelGeometry {
            channel_type: cfd_1d::ChannelType::Straight,
            length: Length::from_base(length),
            cross_section: CrossSection::Circular {
                diameter: Length::from_base(diameter),
            },
            surface: SurfaceProperties {
                roughness: Length::from_base(0.0),
                contact_angle: None,
                surface_energy: None,
                wettability: Wettability::Hydrophilic,
            },
            variations: Vec::new(),
        }),
        resistance_update_policy: ResistanceUpdatePolicy::FlowDependent,
        properties: HashMap::new(),
    };
    network.add_edge_properties(e1, make_props("e1", inlet_diameter, inlet_length));
    network.add_edge_properties(e2, make_props("e2", inlet_diameter, inlet_length));
    network.add_edge_properties(e3, make_props("e3", outlet_diameter, outlet_length));

    let hematocrit_events = vec![
        InletHematocritEvent {
            time: time(0.0),
            node_index: in1.index(),
            hematocrit: fraction(0.45),
        },
        InletHematocritEvent {
            time: time(0.0),
            node_index: in2.index(),
            hematocrit: fraction(0.10),
        },
    ];
    let pressure_events = vec![
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in1.index(),
            pressure: pressure(2.5),
        },
        PressureBoundaryEvent {
            time: time(0.0),
            node_index: in2.index(),
            pressure: pressure(2.0),
        },
    ];

    let states = TransientCompositionSimulator::
        simulate_blood_hematocrit_with_coupled_pressure_events_and_segmented_edge_transport(
            &network,
            hematocrit_events,
        timepoints(&[0.0, 1.0]),
            pressure_events,
            segmented_transport_config(),
        )
        .expect("simulate coupled segmented blood transport");

    let t0_hct = states[0]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);
    let t1_hct = states[1]
        .edge_hematocrit(e3.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);
    let t1_q1 = states[1]
        .edge_flow_rates
        .get(&e1.index())
        .copied()
        .map(|q| q.into_base())
        .expect("branch 1 flow");
    let t1_q2 = states[1]
        .edge_flow_rates
        .get(&e2.index())
        .copied()
        .map(|q| q.into_base())
        .expect("branch 2 flow");
    let t1_h1 = states[1]
        .edge_hematocrit(e1.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);
    let t1_h2 = states[1]
        .edge_hematocrit(e2.index())
        .map(Dimensionless::into_base)
        .unwrap_or(-1.0);
    assert!(t0_hct.abs() < 1.0e-12);
    assert!(t1_hct.is_finite());
    assert!(t1_hct.abs() < 1.0e-12);
    assert!((0.0..=0.45).contains(&t1_h1));
    assert!((0.0..=0.10).contains(&t1_h2));
    assert!(t1_h1 > t1_h2);
    assert!(t1_q1.is_finite() && t1_q1 > 0.0);
    assert!(t1_q2.is_finite() && t1_q2 > 0.0);
}
