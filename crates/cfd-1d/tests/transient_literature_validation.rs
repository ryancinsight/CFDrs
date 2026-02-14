//! Transient literature-validation tests for composition and droplet pipelines.
//!
//! References:
//! - White, F. M. (2015). *Fluid Mechanics* (8th ed.), Eq. for laminar linear
//!   pressure-flow scaling (Hagen-Poiseuille regime): $Q \propto \Delta P$.
//! - Bird, Stewart, Lightfoot (2002). *Transport Phenomena*, continuity/mixing:
//!   ideal instantaneous junction concentration is flow-weighted by inlet streams.
//! - Stone, Stroock, Ajdari (2004). *Engineering flows in small devices*:
//!   particle/droplet advection follows local fluid velocity, $\dot{x}=u=Q/A$.
//!
//! These tests validate transient implementations against analytical expectations.

use cfd_1d::network::{Network, NetworkBuilder};
use cfd_1d::solver::{
    DropletInjection, EdgeFlowEvent, InletCompositionEvent, MixtureComposition,
    PressureBoundaryEvent, TransientCompositionSimulator, TransientDropletSimulator,
};
use std::collections::HashMap;

fn pure(fluid_id: i32) -> MixtureComposition<f64> {
    let mut map = HashMap::new();
    map.insert(fluid_id, 1.0);
    MixtureComposition::new(map)
}

#[test]
fn pressure_event_flow_ratio_matches_hagen_poiseuille_scaling() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_pressure(outlet, 0.0);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: 0.0,
            node_index: inlet.index(),
            pressure: 1000.0,
        },
        PressureBoundaryEvent {
            time: 1.0,
            node_index: inlet.index(),
            pressure: 2000.0,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        vec![0.0, 1.0],
        pressure_events,
    )
    .expect("transient composition");

    let q0 = states[0].edge_flow_rates[&edge.index()].abs();
    let q1 = states[1].edge_flow_rates[&edge.index()].abs();
    assert!((q1 / q0 - 2.0).abs() < 1e-9);
}

#[test]
fn junction_mixing_matches_flow_weighted_mass_conservation() {
    let mut builder = NetworkBuilder::<f64>::new();
    let in1 = builder.add_inlet("in1".to_string());
    let in2 = builder.add_inlet("in2".to_string());
    let j = builder.add_junction("j".to_string());
    let out = builder.add_outlet("out".to_string());

    let _e1 = builder.connect_with_pipe(in1, j, "e1".to_string());
    let _e2 = builder.connect_with_pipe(in2, j, "e2".to_string());
    let e3 = builder.connect_with_pipe(j, out, "e3".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid);
    network.set_pressure(out, 0.0);

    let composition_events = vec![
        InletCompositionEvent {
            time: 0.0,
            node_index: in1.index(),
            mixture: pure(10),
        },
        InletCompositionEvent {
            time: 0.0,
            node_index: in2.index(),
            mixture: pure(20),
        },
    ];

    // With equal branch resistances and p_out=0, choose pressures that keep both inlet
    // streams entering junction and produce analytical 4:1 flow split.
    let pressure_events = vec![
        PressureBoundaryEvent {
            time: 0.0,
            node_index: in1.index(),
            pressure: 3.0,
        },
        PressureBoundaryEvent {
            time: 0.0,
            node_index: in2.index(),
            pressure: 2.0,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        vec![0.0],
        pressure_events,
    )
    .expect("transient composition");

    let mix = states[0]
        .average_fluid_concentrations_in_edge(e3.index())
        .expect("outlet edge composition");

    assert!((mix.get(&10).copied().unwrap_or(0.0) - 0.8).abs() < 1e-9);
    assert!((mix.get(&20).copied().unwrap_or(0.0) - 0.2).abs() < 1e-9);
}

#[test]
fn droplet_advection_matches_velocity_relation_dx_equals_q_over_a_dt() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let network = Network::new(graph, fluid);

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];

    let flow_events = vec![EdgeFlowEvent {
        time: 0.0,
        edge_index: edge.index(),
        flow_rate: 0.5,
    }];

    let injections = vec![DropletInjection {
        droplet_id: 99,
        fluid_id: 1,
        volume: 1e-6,
        injection_time: 0.0,
        channel_index: edge.index(),
        relative_position: 0.1,
    }];

    let states = TransientDropletSimulator::simulate_with_flow_events(
        &network,
        injections,
        composition_events,
        vec![0.0, 0.2],
        flow_events,
    )
    .expect("droplet simulation");

    let pos = states[1].droplets[&99]
        .position
        .as_ref()
        .expect("droplet position")
        .relative_position;

    // area defaults to 1 and length defaults to 1 for this network setup:
    // dx = (Q/A)/L * dt = 0.5 * 0.2 = 0.1, so x = 0.1 + 0.1 = 0.2.
    assert!((pos - 0.2).abs() < 1e-9);
}

#[test]
fn pressure_event_validation_operates_in_laminar_reynolds_range() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let mut network = Network::new(graph, fluid.clone());
    network.set_pressure(outlet, 0.0);

    if let Some(edge_data) = network.graph.edge_weight_mut(edge) {
        edge_data.resistance = 1.0e12;
    }

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];

    let pressure_events = vec![
        PressureBoundaryEvent {
            time: 0.0,
            node_index: inlet.index(),
            pressure: 1000.0,
        },
        PressureBoundaryEvent {
            time: 1.0,
            node_index: inlet.index(),
            pressure: 2000.0,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_pressure_events(
        &network,
        composition_events,
        vec![0.0, 1.0],
        pressure_events,
    )
    .expect("transient composition");

    let diameter: f64 = 100e-6;
    let area = std::f64::consts::PI * (diameter / 2.0_f64).powi(2);

    for state in &states {
        let q = state.edge_flow_rates[&edge.index()].abs();
        let velocity = q / area;
        let re = fluid.density * velocity * diameter / fluid.viscosity;
        assert!(re < 2300.0, "Re={} should remain laminar", re);
        assert!(re < 100.0, "Re={} should remain microfluidic-low", re);
    }
}

#[test]
fn flow_event_validation_operates_in_laminar_reynolds_range() {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("in".to_string());
    let outlet = builder.add_outlet("out".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "e0".to_string());
    let graph = builder.build().expect("graph");

    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
    let network = Network::new(graph, fluid.clone());

    let composition_events = vec![InletCompositionEvent {
        time: 0.0,
        node_index: inlet.index(),
        mixture: pure(1),
    }];

    let flow_events = vec![
        EdgeFlowEvent {
            time: 0.0,
            edge_index: edge.index(),
            flow_rate: 1.0e-9,
        },
        EdgeFlowEvent {
            time: 1.0,
            edge_index: edge.index(),
            flow_rate: 2.0e-9,
        },
    ];

    let states = TransientCompositionSimulator::simulate_with_flow_events(
        &network,
        composition_events,
        vec![0.0, 1.0],
        flow_events,
    )
    .expect("transient composition");

    let diameter: f64 = 100e-6;
    let area = std::f64::consts::PI * (diameter / 2.0_f64).powi(2);

    for state in &states {
        let q = state.edge_flow_rates[&edge.index()].abs();
        let velocity = q / area;
        let re = fluid.density * velocity * diameter / fluid.viscosity;
        assert!(re < 2300.0, "Re={} should remain laminar", re);
        assert!(re < 100.0, "Re={} should remain microfluidic-low", re);
    }
}
