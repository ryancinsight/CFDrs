//! Comprehensive turbulence model performance benchmarks
//!
//! Benchmarks CPU vs GPU performance across turbulence model hierarchy:
//! - RANS models: k-ε, k-ω SST, Spalart-Allmaras
//! - LES models: Smagorinsky LES, DES
//! - Performance scaling with grid size and model complexity

use cfd_2d::physics::turbulence::*;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

/// Benchmark RANS model performance (k-ε, k-ω SST, SA)
fn bench_rans_models(c: &mut Criterion) {
    let mut group = c.benchmark_group("rans_models");

    // Grid sizes to test scaling
    let grid_sizes = vec![(32, 32), (64, 64), (128, 128)];

    for (nx, ny) in grid_sizes {
        // k-ε model benchmark
        group.bench_function(format!("k_epsilon_{}x{}", nx, ny), |b| {
            let mut k_model = KEpsilonModel::<f64>::new(nx, ny);
            let mut k_field = vec![0.1; nx * ny];
            let mut eps_field = vec![0.01; nx * ny];
            let velocity_field = vec![nalgebra::Vector2::new(1.0, 0.0); nx * ny];

            b.iter(|| {
                k_model
                    .update(
                        &mut k_field,
                        &mut eps_field,
                        &velocity_field,
                        black_box(1.0),    // density
                        black_box(1.5e-5), // viscosity
                        black_box(1e-4),   // dt
                        black_box(0.01),   // dx
                        black_box(0.01),   // dy
                    )
                    .unwrap();
            });
        });

        // k-ω SST model benchmark
        group.bench_function(format!("k_omega_sst_{}x{}", nx, ny), |b| {
            let mut k_model = KOmegaSSTModel::<f64>::new(nx, ny);
            let mut k_field = vec![0.1; nx * ny];
            let mut omega_field = vec![10.0; nx * ny];
            let velocity_field = vec![nalgebra::Vector2::new(1.0, 0.0); nx * ny];

            b.iter(|| {
                k_model
                    .update(
                        &mut k_field,
                        &mut omega_field,
                        &velocity_field,
                        black_box(1.0),
                        black_box(1.5e-5),
                        black_box(1e-4),
                        black_box(0.01),
                        black_box(0.01),
                    )
                    .unwrap();
            });
        });

        // Spalart-Allmaras model benchmark
        group.bench_function(format!("spalart_allmaras_{}x{}", nx, ny), |b| {
            let sa_model = SpalartAllmaras::<f64>::new(nx, ny);
            let mut nu_tilde_field = vec![1e-4; nx * ny];
            let velocity_field = vec![nalgebra::Vector2::new(1.0, 0.0); nx * ny];

            b.iter(|| {
                sa_model
                    .update(
                        &mut nu_tilde_field,
                        &velocity_field,
                        black_box(1.5e-5),
                        black_box(1e-4),
                        black_box(0.01),
                        black_box(0.01),
                    )
                    .unwrap();
            });
        });
    }

    group.finish();
}

/// Benchmark LES model performance (Smagorinsky LES, DES)
fn bench_les_models(c: &mut Criterion) {
    let mut group = c.benchmark_group("les_models");

    let grid_sizes = vec![(32, 32), (64, 64)];

    for (nx, ny) in grid_sizes {
        // Smagorinsky LES benchmark
        group.bench_function(format!("smagorinsky_les_{}x{}", nx, ny), |b| {
            let config = les_smagorinsky::SmagorinskyConfig {
                smagorinsky_constant: 0.1,
                dynamic_procedure: false,
                wall_damping: false,
                van_driest_constant: 0.0,
                min_sgs_viscosity: 1e-10,
                use_gpu: false, // CPU benchmark
            };
            let mut les_model = SmagorinskyLES::new(nx, ny, 0.01, 0.01, config);
            let velocity_u = nalgebra::DMatrix::from_element(nx, ny, 1.0);
            let velocity_v = nalgebra::DMatrix::from_element(nx, ny, 0.0);
            let pressure = nalgebra::DMatrix::zeros(nx, ny);

            b.iter(|| {
                les_model
                    .update(
                        &velocity_u,
                        &velocity_v,
                        &pressure,
                        black_box(1.0),
                        black_box(1.5e-5),
                        black_box(1e-4),
                        black_box(0.01),
                        black_box(0.01),
                    )
                    .unwrap();
            });
        });

        // DES benchmark
        group.bench_function(format!("des_{}x{}", nx, ny), |b| {
            let config = des::DESConfig {
                variant: des::DESVariant::DDES,
                des_constant: 0.65,
                max_sgs_ratio: 0.5,
                rans_viscosity: 1.5e-5,
                use_gpu: false, // CPU benchmark
            };
            let mut des_model = DetachedEddySimulation::new(nx, ny, 0.01, 0.01, config, &[]);
            let velocity_u = nalgebra::DMatrix::from_element(nx, ny, 1.0);
            let velocity_v = nalgebra::DMatrix::from_element(nx, ny, 0.0);
            let pressure = nalgebra::DMatrix::zeros(nx, ny);

            b.iter(|| {
                des_model
                    .update(
                        &velocity_u,
                        &velocity_v,
                        &pressure,
                        black_box(1.0),
                        black_box(1.5e-5),
                        black_box(1e-4),
                        black_box(0.01),
                        black_box(0.01),
                    )
                    .unwrap();
            });
        });
    }

    group.finish();
}

/// Benchmark turbulence model initialization overhead
fn bench_model_initialization(c: &mut Criterion) {
    let mut group = c.benchmark_group("model_initialization");

    group.bench_function("k_epsilon_init", |b| {
        b.iter(|| {
            black_box(KEpsilonModel::<f64>::new(64, 64));
        });
    });

    group.bench_function("k_omega_sst_init", |b| {
        b.iter(|| {
            black_box(KOmegaSSTModel::<f64>::new(64, 64));
        });
    });

    group.bench_function("spalart_allmaras_init", |b| {
        b.iter(|| {
            black_box(SpalartAllmaras::<f64>::new(64, 64));
        });
    });

    group.bench_function("smagorinsky_les_init", |b| {
        let config = les_smagorinsky::SmagorinskyConfig {
            smagorinsky_constant: 0.1,
            dynamic_procedure: false,
            wall_damping: false,
            van_driest_constant: 0.0,
            min_sgs_viscosity: 1e-10,
            use_gpu: false,
        };
        b.iter(|| {
            black_box(SmagorinskyLES::new(64, 64, 0.01, 0.01, config.clone()));
        });
    });

    group.bench_function("des_init", |b| {
        let config = des::DESConfig {
            variant: des::DESVariant::DDES,
            des_constant: 0.65,
            max_sgs_ratio: 0.5,
            rans_viscosity: 1.5e-5,
            use_gpu: false,
        };
        b.iter(|| {
            black_box(DetachedEddySimulation::new(
                64,
                64,
                0.01,
                0.01,
                config.clone(),
                &[],
            ));
        });
    });

    group.finish();
}

/// Benchmark turbulence property calculations
fn bench_turbulence_properties(c: &mut Criterion) {
    let mut group = c.benchmark_group("turbulence_properties");

    // Turbulent viscosity calculations
    group.bench_function("k_epsilon_viscosity", |b| {
        let model = KEpsilonModel::<f64>::new(1, 1);
        b.iter(|| {
            black_box(model.turbulent_viscosity(0.1, 0.01, 1.0));
        });
    });

    group.bench_function("k_omega_sst_viscosity", |b| {
        let model = KOmegaSSTModel::<f64>::new(1, 1);
        b.iter(|| {
            black_box(model.turbulent_viscosity_with_limiter(0.1, 10.0, 1.0, 100.0, 1.0));
        });
    });

    group.bench_function("spalart_allmaras_viscosity", |b| {
        let model = SpalartAllmaras::<f64>::new(1, 1);
        b.iter(|| {
            black_box(model.eddy_viscosity(1e-4, 1.5e-5));
        });
    });

    group.finish();
}

/// Benchmark memory usage patterns
fn bench_memory_usage(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_usage");

    group.bench_function("field_allocation_64x64", |b| {
        b.iter(|| {
            let _k = vec![0.1f64; 64 * 64];
            let _eps = vec![0.01f64; 64 * 64];
            let _velocity = vec![nalgebra::Vector2::<f64>::zeros(); 64 * 64];
            black_box(());
        });
    });

    group.bench_function("matrix_allocation_64x64", |b| {
        b.iter(|| {
            let _velocity_u = nalgebra::DMatrix::<f64>::zeros(64, 64);
            let _velocity_v = nalgebra::DMatrix::<f64>::zeros(64, 64);
            let _pressure = nalgebra::DMatrix::<f64>::zeros(64, 64);
            black_box(());
        });
    });

    group.finish();
}

criterion_group! {
    name = turbulence_benchmarks;
    config = Criterion::default()
        .sample_size(10)
        .measurement_time(std::time::Duration::from_secs(5));
    targets = bench_rans_models, bench_les_models, bench_model_initialization,
             bench_turbulence_properties, bench_memory_usage
}

criterion_main!(turbulence_benchmarks);
