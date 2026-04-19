use cfd_2d::{
    grid::StructuredGrid2D,
    network::{solve_reference_trace, Network2dBuilderSink},
    solvers::lbm::{LbmConfig, LbmSolver},
};
use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
use cfd_schematics::interface::presets::{primitive_selective_split_tree_rect, venturi_rect};
use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use nalgebra::Vector2;
use std::collections::HashMap;

fn benchmark_lbm_solver(c: &mut Criterion) {
    let mut group = c.benchmark_group("lbm_solver");

    for size in [32, 64, 128].iter() {
        let grid = StructuredGrid2D::new(*size, *size, 0.0, 1.0, 0.0, 1.0).unwrap();
        let config = LbmConfig::default();
        let mut solver = LbmSolver::new(config, &grid);

        // Initialize with some density and velocity
        solver
            .initialize(
                |_x: f64, _y: f64| 1.0,
                |_x: f64, _y: f64| Vector2::new(0.0, 0.0),
            )
            .unwrap();

        group.bench_with_input(BenchmarkId::new("step", size), size, |b, _| {
            b.iter(|| {
                let boundaries = HashMap::new();
                let _: () = solver.step(&boundaries).unwrap();
                black_box(())
            })
        });
    }

    group.finish();
}

fn benchmark_grid_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("grid_creation");

    for size in [50, 100, 200].iter() {
        group.bench_with_input(
            BenchmarkId::new("structured_grid", size),
            size,
            |b, &size| {
                b.iter(|| black_box(StructuredGrid2D::new(size, size, 0.0, 1.0, 0.0, 1.0).unwrap()))
            },
        );
    }

    group.finish();
}

fn with_provenance(
    mut blueprint: cfd_schematics::domain::model::NetworkBlueprint,
) -> cfd_schematics::domain::model::NetworkBlueprint {
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(
            cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
        );
    blueprint
}

fn benchmark_network_pipeline(c: &mut Criterion) {
    let mut group = c.benchmark_group("network_pipeline");

    let cases = [
        (
            "venturi",
            with_provenance(venturi_rect(
                "bench_venturi",
                2.0e-3,
                0.7e-3,
                1.0e-3,
                1.8e-3,
            )),
            24usize,
            10usize,
            1.0e-6_f64,
        ),
        (
            "selective",
            with_provenance(primitive_selective_split_tree_rect(
                "bench_selective",
                (127.76, 85.47),
                &[PrimitiveSelectiveSplitKind::Tri],
                2.0e-3,
                0.45,
                0.45,
                0.68,
                0.5e-3,
                1.0e-3,
                1.0e-3,
                false,
                0,
                None,
            )),
            24usize,
            10usize,
            1.0e-6_f64,
        ),
    ];

    for (label, blueprint, grid_nx, grid_ny, q_total) in cases {
        group.bench_with_input(
            BenchmarkId::new("reference_trace", label),
            &blueprint,
            |b, blueprint| {
                b.iter(|| {
                    black_box(
                        solve_reference_trace::<f64>(blueprint, 1060.0, 3.5e-3, q_total)
                            .expect("reference trace"),
                    )
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("build_with_projection", label),
            &blueprint,
            |b, blueprint| {
                b.iter(|| {
                    let blood = BloodModel::Newtonian(3.5e-3_f64);
                    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, grid_nx, grid_ny);
                    black_box(sink.build(blueprint).expect("network build"))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("core_solve", label),
            &blueprint,
            |b, blueprint| {
                b.iter_batched(
                    || {
                        let blood = BloodModel::Newtonian(3.5e-3_f64);
                        let sink =
                            Network2dBuilderSink::new(blood, 1060.0, q_total, grid_nx, grid_ny);
                        sink.build(blueprint).expect("network build")
                    },
                    |mut network| black_box(network.solve_all(1e-6).expect("solve all")),
                    BatchSize::SmallInput,
                )
            },
        );

        group.bench_with_input(
            BenchmarkId::new("projected_solve", label),
            &blueprint,
            |b, blueprint| {
                b.iter_batched(
                    || {
                        let blood = BloodModel::Newtonian(3.5e-3_f64);
                        let sink =
                            Network2dBuilderSink::new(blood, 1060.0, q_total, grid_nx, grid_ny);
                        sink.build(blueprint).expect("network build")
                    },
                    |mut network| {
                        black_box(network.solve_projected(1e-6).expect("projected solve"))
                    },
                    BatchSize::SmallInput,
                )
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    benchmark_lbm_solver,
    benchmark_grid_creation,
    benchmark_network_pipeline
);
criterion_main!(benches);
