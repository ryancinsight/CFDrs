//! Criterion benchmarks for convergence monitoring and analysis
//!
//! These benchmarks measure performance of convergence detection algorithms.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use cfd_validation::convergence::{ConvergenceMonitor, GridConvergenceIndex};

fn benchmark_convergence_monitor(c: &mut Criterion) {
    let mut group = c.benchmark_group("convergence_monitor");
    
    for size in [10, 100, 1000].iter() {
        group.throughput(Throughput::Elements(*size as u64));
        
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, size);
                let mut error = 1.0;
                
                for _ in 0..size {
                    monitor.update(black_box(error));
                    error *= 0.9;
                }
                
                black_box(monitor.check_status())
            });
        });
    }
    
    group.finish();
}

fn benchmark_gci_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("gci_calculation");
    
    let gci = GridConvergenceIndex::<f64>::new(3, 2.0, 2.0);
    let values = vec![
        (1.0, 1.01),
        (10.0, 10.1),
        (100.0, 100.5),
        (1000.0, 1005.0),
    ];
    
    for (f_fine, f_coarse) in values {
        group.bench_with_input(
            BenchmarkId::new("compute_fine", f_fine),
            &(f_fine, f_coarse),
            |b, &(f_fine, f_coarse)| {
                b.iter(|| {
                    black_box(gci.compute_fine(black_box(f_fine), black_box(f_coarse)))
                });
            },
        );
    }
    
    group.finish();
}

fn benchmark_convergence_status_check(c: &mut Criterion) {
    let mut group = c.benchmark_group("convergence_status_check");
    
    // Pre-populate monitor with history
    for history_size in [10, 50, 100].iter() {
        group.bench_with_input(
            BenchmarkId::from_parameter(history_size),
            history_size,
            |b, &size| {
                let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 1000);
                let mut error = 1.0;
                
                for _ in 0..size {
                    monitor.update(error);
                    error *= 0.95;
                }
                
                b.iter(|| {
                    black_box(monitor.check_status())
                });
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    benchmark_convergence_monitor,
    benchmark_gci_calculation,
    benchmark_convergence_status_check
);
criterion_main!(benches);
