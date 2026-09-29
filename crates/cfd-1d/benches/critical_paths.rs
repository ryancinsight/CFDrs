#![allow(missing_docs)]
//! Critical-path enumeration cost of `ResistanceAnalyzer` on a ladder network.
//!
//! A ladder of `stages` junction hops with `parallel` pipes per hop has one
//! node path and `parallel^stages` edge paths, so the benchmark isolates the
//! per-edge-path work of the analyzer rather than graph path search.
use aequitas::systems::si::quantities::{Area, HydraulicResistance, Length};
use cfd_1d::solver::analysis::analyzers::{NetworkAnalyzer, ResistanceAnalyzer};
use cfd_1d::{ComponentType, EdgeProperties, Network, NetworkBuilder, ResistanceUpdatePolicy};
use cfd_core::physics::fluid::database::water_20c;
use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use std::collections::HashMap;

fn ladder_network(stages: usize, parallel: usize) -> Network<f64> {
    let mut builder = NetworkBuilder::new();
    let mut previous = builder.add_inlet("inlet".to_string());
    let mut pipes = Vec::with_capacity(stages * parallel);
    for stage in 0..stages {
        let next = if stage + 1 == stages {
            builder.add_outlet("outlet".to_string())
        } else {
            builder.add_junction(format!("j{stage}"))
        };
        for branch in 0..parallel {
            let id = format!("s{stage}b{branch}");
            pipes.push((
                builder.connect_with_pipe(previous, next, id.clone()),
                id,
                branch,
            ));
        }
        previous = next;
    }
    let graph = builder.build().expect("ladder network must validate");
    let mut network = Network::new(graph, water_20c::<f64>().expect("water properties exist"));
    for (edge, id, branch) in pipes {
        // Distinct diameters per branch give distinct Hagen-Poiseuille resistances.
        let diameter = 1.0e-3 * (1.0 + 0.25 * branch as f64);
        network.add_edge_properties(
            edge,
            EdgeProperties {
                id,
                component_type: ComponentType::Pipe,
                length: Length::from_base(0.01),
                area: Area::from_base(std::f64::consts::FRAC_PI_4 * diameter * diameter),
                hydraulic_diameter: Some(Length::from_base(diameter)),
                resistance: HydraulicResistance::from_base(1.0),
                geometry: None,
                resistance_update_policy: ResistanceUpdatePolicy::FlowInvariant,
                properties: HashMap::new(),
            },
        );
    }
    network
}

fn bench_critical_paths(c: &mut Criterion) {
    let mut group = c.benchmark_group("resistance_critical_paths");
    for (stages, parallel) in [(5, 2), (3, 3)] {
        let network = ladder_network(stages, parallel);
        group.bench_with_input(
            BenchmarkId::from_parameter(format!("{stages}x{parallel}")),
            &network,
            |b, network| {
                let mut analyzer = ResistanceAnalyzer::<f64>::new();
                b.iter(|| {
                    analyzer
                        .analyze(black_box(network))
                        .expect("ladder analysis must succeed")
                });
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_critical_paths);
criterion_main!(benches);
