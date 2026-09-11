#![cfg_attr(test, expect(clippy::print_stdout, reason = "test/validation output"))]
//! Explicit opt-in allocation-instrumentation harness for the tracking tests.
//!
//! `TrackingAllocator` is a process-global instrument: its counts are only
//! meaningful while it is the sole `#[global_allocator]` of the process
//! running the measurement. These tests therefore live in their own
//! integration-test binary, which installs the allocator explicitly — the
//! CFDRS-GA-001 contract ("allocator only in explicitly opted-in bench/bin
//! targets"). No other test, bench, example, or downstream binary in the
//! crate's link graph inherits instrumentation it did not request.

use cfd_validation::benchmarking::{CfdMemoryProfiler, MemoryProfiler, TrackingAllocator};

#[global_allocator]
static TRACKING_ALLOCATOR: TrackingAllocator = TrackingAllocator::new();

#[test]
fn profiler_reads_counts_from_the_installed_tracking_allocator() {
    let profiler = MemoryProfiler::new(TRACKING_ALLOCATOR.stats());

    let (result, stats) = profiler
        .profile_closure(|| {
            let data = vec![1.0f64; 1000];
            data.iter().sum::<f64>()
        })
        .expect("invariant: the test profiling session starts before it stops");

    assert!(result > 0.0);
    assert!(stats.total_allocated > 0);
    assert!(stats.peak_allocated >= stats.current_allocated);
}

#[test]
fn cfd_memory_suite_records_allocations_from_the_installed_allocator() {
    let cfd_profiler = CfdMemoryProfiler::new(TRACKING_ALLOCATOR.stats());
    let results = cfd_profiler
        .run_memory_suite()
        .expect("invariant: the test profiling suite has valid inputs");

    assert!(!results.is_empty());
    for (name, stats) in results {
        println!("{name}: {stats}");
        assert!(stats.total_allocated > 0);
    }
}
