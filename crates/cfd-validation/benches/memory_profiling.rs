#![allow(missing_docs)]
//! Explicit allocation-instrumentation benchmark for CFD validation workloads.

use cfd_validation::benchmarking::{CfdMemoryProfiler, TrackingAllocator};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

#[global_allocator]
static GLOBAL_ALLOCATOR: TrackingAllocator = TrackingAllocator::new();

fn benchmark_memory_suite(c: &mut Criterion) {
    c.bench_function("cfd_memory_suite", |b| {
        b.iter(|| {
            let profiler = CfdMemoryProfiler::new(GLOBAL_ALLOCATOR.stats());
            let results = profiler
                .run_memory_suite()
                .expect("benchmark invariant: CFD memory suite inputs are valid");
            black_box(results);
        });
    });
}

criterion_group!(benches, benchmark_memory_suite);
criterion_main!(benches);
